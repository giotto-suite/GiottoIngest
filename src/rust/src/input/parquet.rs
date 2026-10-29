//! A parquet source read the way GiottoDisk's lazy query reads it: row
//! groups in file order, the reader's filters applied, coordinates as
//! doubles with y optionally negated.

use std::fs::File;

use std::sync::Arc;

use arrow::array::{new_null_array, Array, ArrayRef, BooleanArray, Float64Array};
use arrow::compute::kernels::cmp;
use arrow::compute::{and, cast, filter_record_batch, or};
use arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ProjectionMask;
use rayon::prelude::*;

pub struct Filter {
    pub col: String,
    pub op: String,
    /// one value, or the set for `in`
    pub values: Vec<f64>,
}

pub struct Spec {
    pub files: Vec<String>,
    pub x: String,
    pub y: String,
    pub flip_y: bool,
    pub filters: Vec<Filter>,
}

/// One unit of parallel decoding: a row group of a file.
struct Unit {
    file: usize,
    rg: usize,
}

fn schema_of(path: &str) -> Result<SchemaRef, String> {
    let f = File::open(path).map_err(|e| format!("{path}: {e}"))?;
    let b = ParquetRecordBatchReaderBuilder::try_new(f).map_err(|e| format!("{path}: {e}"))?;
    Ok(b.schema().clone())
}

pub fn source_schema(spec: &Spec) -> Result<SchemaRef, String> {
    schema_of(spec.files.first().ok_or("no input files")?)
}

/// Whether the x and y columns are float32 (comparisons then happen in f32).
pub fn xy_f32(spec: &Spec) -> Result<(bool, bool), String> {
    let s = source_schema(spec)?;
    let is32 = |c: &str| -> Result<bool, String> {
        Ok(s.field_with_name(c).map_err(|_| format!("no column `{c}`"))?.data_type() == &DataType::Float32)
    };
    Ok((is32(&spec.x)?, is32(&spec.y)?))
}

fn units(spec: &Spec) -> Result<Vec<Unit>, String> {
    let mut out = Vec::new();
    for (fi, path) in spec.files.iter().enumerate() {
        let f = File::open(path).map_err(|e| format!("{path}: {e}"))?;
        let b = ParquetRecordBatchReaderBuilder::try_new(f).map_err(|e| format!("{path}: {e}"))?;
        for rg in 0..b.metadata().num_row_groups() {
            out.push(Unit { file: fi, rg });
        }
    }
    Ok(out)
}

/// Reads `cols` of one row group as the first file's schema has them, the
/// way R's `open_dataset()` reads several files: later files' columns are
/// cast to the first file's types, and columns a later file lacks are null.
fn read_unit(spec: &Spec, u: &Unit, cols: &[String], target: &Schema) -> Result<Vec<RecordBatch>, String> {
    let path = &spec.files[u.file];
    let f = File::open(path).map_err(|e| format!("{path}: {e}"))?;
    let b = ParquetRecordBatchReaderBuilder::try_new(f).map_err(|e| format!("{path}: {e}"))?;
    let sch = b.schema().clone();
    let mut roots: Vec<usize> = cols.iter().filter_map(|c| sch.index_of(c).ok()).collect();
    roots.sort_unstable();
    let mask = ProjectionMask::roots(b.parquet_schema(), roots);
    let rdr = b
        .with_row_groups(vec![u.rg])
        .with_projection(mask)
        .with_batch_size(1 << 16)
        .build()
        .map_err(|e| format!("{path}: {e}"))?;
    let out_fields: Vec<Field> = cols.iter()
        .map(|c| target.field_with_name(c).cloned().map_err(|_| format!("no column `{c}`")))
        .collect::<Result<_, _>>()?;
    let out_schema = Arc::new(Schema::new(out_fields));
    rdr.map(|r| {
        let b = r.map_err(|e| format!("{path}: {e}"))?;
        let n = b.num_rows();
        let arrays: Vec<ArrayRef> = out_schema.fields().iter().map(|f| match b.column_by_name(f.name()) {
            Some(a) if a.data_type() == f.data_type() => Ok(a.clone()),
            Some(a) => cast(a, f.data_type()).map_err(|e| format!("{path}: cast `{}`: {e}", f.name())),
            None => Ok(new_null_array(f.data_type(), n)),
        }).collect::<Result<_, _>>()?;
        RecordBatch::try_new(out_schema.clone(), arrays).map_err(|e| format!("{path}: {e}"))
    }).collect()
}

fn as_f64(a: &ArrayRef) -> Result<Float64Array, String> {
    let c = cast(a, &DataType::Float64).map_err(|e| e.to_string())?;
    Ok(c.as_any().downcast_ref::<Float64Array>().unwrap().clone())
}

/// The reader's filters, as R's arrow query evaluates them: numeric
/// comparison in double, NA rows dropped.
fn apply_filters(spec: &Spec, b: RecordBatch) -> Result<RecordBatch, String> {
    if spec.filters.is_empty() {
        return Ok(b);
    }
    let mut mask: Option<BooleanArray> = None;
    for f in &spec.filters {
        let col = b.column_by_name(&f.col).ok_or(format!("no column `{}`", f.col))?;
        // numeric only: arrow would otherwise parse strings as numbers
        if !(col.data_type().is_integer() || col.data_type().is_floating()) {
            return Err(format!("unsupported: filter column `{}` is {}, not numeric",
                f.col, col.data_type()));
        }
        let v_col = as_f64(col)?;
        // arrow casts the scalar to a float32 column's type before comparing
        let m = if f.op == "in" {
            // R's %in% on a numeric column: equality with any value, NA dropped
            let mut m: Option<BooleanArray> = None;
            for v in &f.values {
                let val = if col.data_type() == &DataType::Float32 { *v as f32 as f64 } else { *v };
                let e = cmp::eq(&v_col, &Float64Array::new_scalar(val)).map_err(|e| e.to_string())?;
                m = Some(match m { None => e, Some(p) => or(&p, &e).map_err(|e| e.to_string())? });
            }
            m.unwrap_or_else(|| BooleanArray::from(vec![false; v_col.len()]))
        } else {
            let one = *f.values.first().ok_or(format!("filter on `{}` has no value", f.col))?;
            let val = if col.data_type() == &DataType::Float32 { one as f32 as f64 } else { one };
            let s = Float64Array::new_scalar(val);
            match f.op.as_str() {
                ">=" => cmp::gt_eq(&v_col, &s),
                ">" => cmp::gt(&v_col, &s),
                "<=" => cmp::lt_eq(&v_col, &s),
                "<" => cmp::lt(&v_col, &s),
                "==" => cmp::eq(&v_col, &s),
                "!=" => cmp::neq(&v_col, &s),
                o => return Err(format!("unsupported filter op `{o}`")),
            }
            .map_err(|e| e.to_string())?
        };
        mask = Some(match mask {
            None => m,
            Some(p) => and(&p, &m).map_err(|e| e.to_string())?,
        });
    }
    filter_record_batch(&b, &mask.unwrap()).map_err(|e| e.to_string())
}

/// Coordinates of a (filtered) batch as doubles, y negated when flipping.
/// A null coordinate is NaN, so it is in no tile and not in the extent, as
/// R's comparisons and `min(na.rm = TRUE)` treat NA; the row still counts in
/// the source's `n`, as in R's `nrow()`.
pub fn xy(spec: &Spec, b: &RecordBatch) -> Result<(Vec<f64>, Vec<f64>), String> {
    let x = as_f64(b.column_by_name(&spec.x).ok_or(format!("no column `{}`", spec.x))?)?;
    let y = as_f64(b.column_by_name(&spec.y).ok_or(format!("no column `{}`", spec.y))?)?;
    let vals = |a: &Float64Array, neg: bool| -> Vec<f64> {
        let s = if neg { -1.0 } else { 1.0 };
        if a.null_count() == 0 {
            a.values().iter().map(|v| v * s).collect()
        } else {
            a.iter().map(|v| v.map_or(f64::NAN, |v| v * s)).collect()
        }
    };
    Ok((vals(&x, false), vals(&y, spec.flip_y)))
}

/// Calls `f` on every filtered batch, in source order. Row groups are decoded
/// `threads` at a time in parallel; `f` itself runs sequentially.
pub fn for_each_batch(
    spec: &Spec,
    cols: &[String],
    threads: usize,
    mut f: impl FnMut(RecordBatch) -> Result<(), String>,
) -> Result<(), String> {
    let mut need: Vec<String> = cols.to_vec();
    for c in [&spec.x, &spec.y].into_iter().chain(spec.filters.iter().map(|f| &f.col)) {
        if !need.contains(c) {
            need.push(c.clone());
        }
    }
    let us = units(spec)?;
    let target = source_schema(spec)?;
    let window = (threads * 2).max(1);
    for chunk in us.chunks(window) {
        crate::core::cancel::check()?;
        let decoded: Vec<Result<Vec<RecordBatch>, String>> = chunk
            .par_iter()
            .map(|u| {
                read_unit(spec, u, &need, &target)?
                    .into_iter()
                    .map(|b| apply_filters(spec, b))
                    .collect()
            })
            .collect();
        for r in decoded {
            for b in r? {
                crate::core::cancel::check()?;
                if b.num_rows() > 0 {
                    f(b)?;
                }
            }
        }
    }
    Ok(())
}

/// Maps every row group to a partial result in parallel and folds the
/// partials; for order-free passes (extent, counts).
pub fn reduce_units<T: Send>(
    spec: &Spec,
    cols: &[String],
    init: impl Fn() -> T + Sync,
    each: impl Fn(&mut T, RecordBatch) -> Result<(), String> + Sync,
    merge: impl Fn(T, T) -> T + Sync + Send,
) -> Result<T, String> {
    let mut need: Vec<String> = cols.to_vec();
    for c in [&spec.x, &spec.y].into_iter().chain(spec.filters.iter().map(|f| &f.col)) {
        if !need.contains(c) {
            need.push(c.clone());
        }
    }
    let us = units(spec)?;
    let target = source_schema(spec)?;
    us.par_iter()
        .map(|u| -> Result<T, String> {
            crate::core::cancel::check()?;
            let mut acc = init();
            for b in read_unit(spec, u, &need, &target)? {
                let b = apply_filters(spec, b)?;
                if b.num_rows() > 0 {
                    each(&mut acc, b)?;
                }
            }
            Ok(acc)
        })
        .try_reduce(&init, |a, b| Ok(merge(a, b)))
}
