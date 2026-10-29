//! A parquet source read the way GiottoDisk's lazy query reads it: row
//! groups in file order, the reader's filters applied, coordinates as
//! doubles with y optionally negated.

use std::fs::File;

use arrow::array::{ArrayRef, BooleanArray, Float64Array};
use arrow::compute::kernels::cmp;
use arrow::compute::{and, cast, filter_record_batch};
use arrow::datatypes::{DataType, SchemaRef};
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ProjectionMask;
use rayon::prelude::*;

pub struct Filter {
    pub col: String,
    pub op: String,
    pub value: f64,
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

fn read_unit(spec: &Spec, u: &Unit, cols: &[String]) -> Result<Vec<RecordBatch>, String> {
    let path = &spec.files[u.file];
    let f = File::open(path).map_err(|e| format!("{path}: {e}"))?;
    let b = ParquetRecordBatchReaderBuilder::try_new(f).map_err(|e| format!("{path}: {e}"))?;
    let sch = b.schema().clone();
    let roots: Vec<usize> = cols
        .iter()
        .map(|c| sch.index_of(c).map_err(|_| format!("{path}: no column `{c}`")))
        .collect::<Result<_, _>>()?;
    let mask = ProjectionMask::roots(b.parquet_schema(), roots);
    let rdr = b
        .with_row_groups(vec![u.rg])
        .with_projection(mask)
        .with_batch_size(1 << 16)
        .build()
        .map_err(|e| format!("{path}: {e}"))?;
    rdr.map(|r| r.map_err(|e| format!("{path}: {e}"))).collect()
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
        let v = as_f64(col)?;
        // arrow casts the scalar to a float32 column's type before comparing
        let val = if col.data_type() == &DataType::Float32 { f.value as f32 as f64 } else { f.value };
        let s = Float64Array::new_scalar(val);
        let m = match f.op.as_str() {
            ">=" => cmp::gt_eq(&v, &s),
            ">" => cmp::gt(&v, &s),
            "<=" => cmp::lt_eq(&v, &s),
            "<" => cmp::lt(&v, &s),
            "==" => cmp::eq(&v, &s),
            "!=" => cmp::neq(&v, &s),
            o => return Err(format!("unsupported filter op `{o}`")),
        }
        .map_err(|e| e.to_string())?;
        mask = Some(match mask {
            None => m,
            Some(p) => and(&p, &m).map_err(|e| e.to_string())?,
        });
    }
    filter_record_batch(&b, &mask.unwrap()).map_err(|e| e.to_string())
}

/// Coordinates of a (filtered) batch as doubles, y negated when flipping.
pub fn xy(spec: &Spec, b: &RecordBatch) -> Result<(Vec<f64>, Vec<f64>), String> {
    let x = as_f64(b.column_by_name(&spec.x).ok_or(format!("no column `{}`", spec.x))?)?;
    let y = as_f64(b.column_by_name(&spec.y).ok_or(format!("no column `{}`", spec.y))?)?;
    let xs: Vec<f64> = x.values().to_vec();
    let ys: Vec<f64> = if spec.flip_y { y.values().iter().map(|v| -v).collect() } else { y.values().to_vec() };
    Ok((xs, ys))
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
    let window = (threads * 2).max(1);
    for chunk in us.chunks(window) {
        let decoded: Vec<Result<Vec<RecordBatch>, String>> = chunk
            .par_iter()
            .map(|u| {
                read_unit(spec, u, &need)?
                    .into_iter()
                    .map(|b| apply_filters(spec, b))
                    .collect()
            })
            .collect();
        for r in decoded {
            for b in r? {
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
    us.par_iter()
        .map(|u| -> Result<T, String> {
            let mut acc = init();
            for b in read_unit(spec, u, &need)? {
                let b = apply_filters(spec, b)?;
                if b.num_rows() > 0 {
                    each(&mut acc, b)?;
                }
            }
            Ok(acc)
        })
        .try_reduce(&init, |a, b| Ok(merge(a, b)))
}
