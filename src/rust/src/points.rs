//! Point tile store passes: extent, per-depth tile counts, and the routed
//! write. The write reproduces GiottoDisk's per-tile path
//! (`storeWrite(parquetGeomStore, data.frame)` -> `.store_write_points`):
//! rows in source order numbered 1..n, Morton order, WKB points, and the
//! column types an R round trip produces.

use std::collections::HashMap;
use std::fs::{create_dir_all, File};
use std::sync::Arc;

use arrow::array::{
    Array, ArrayRef, AsArray, DictionaryArray, Float64Array, Int16Array, Int32Array, Int8Array,
    RecordBatch, StringArray, UInt32Array,
};
use arrow::compute::{cast, concat_batches, take};
use arrow::datatypes::{DataType, Field, Int16Type, Int32Type, Int8Type, Schema, SchemaRef};
use arrow::ipc::reader::FileReader;
use arrow::ipc::writer::FileWriter;
use rayon::prelude::*;

use crate::core::morton::order;
use crate::core::tile_write::{base_columns, write_parquet};
use crate::input::parquet::{for_each_batch, reduce_units, source_schema, xy, xy_f32, Spec};
use crate::tiles::Tiles;

pub struct Extent {
    pub n: f64,
    pub e: [f64; 4],
}

pub fn extent(spec: &Spec, _threads: usize) -> Result<Extent, String> {
    let init = || Extent { n: 0.0, e: [f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY] };
    reduce_units(spec, &[], init, |acc, b| {
        let (xs, ys) = xy(spec, &b)?;
        acc.n += xs.len() as f64;
        for (x, y) in xs.iter().zip(&ys) {
            if !x.is_nan() { acc.e[0] = acc.e[0].min(*x); acc.e[1] = acc.e[1].max(*x); }
            if !y.is_nan() { acc.e[2] = acc.e[2].min(*y); acc.e[3] = acc.e[3].max(*y); }
        }
        Ok(())
    }, |a, c| Extent { n: a.n + c.n, e: [a.e[0].min(c.e[0]), a.e[1].max(c.e[1]), a.e[2].min(c.e[2]), a.e[3].max(c.e[3])] })
}

pub fn count_tiles(spec: &Spec, bounds: Vec<[f64; 4]>, _threads: usize) -> Result<Vec<f64>, String> {
    let (fx, fy) = xy_f32(spec)?;
    let tiles = Tiles::new(bounds, fx, fy);
    let nt = tiles.len();
    reduce_units(spec, &[], || vec![0f64; nt], |acc, b| {
        let (xs, ys) = xy(spec, &b)?;
        let mut hit = Vec::new();
        for (x, y) in xs.iter().zip(&ys) {
            tiles.matches(*x, *y, &mut hit);
            for k in &hit { acc[*k as usize] += 1.0; }
        }
        Ok(())
    }, |mut a, c| { for (u, v) in a.iter_mut().zip(c) { *u += v; } a })
}

/// How a source column is carried through the spill and typed at write
/// time, following the R round trip GiottoDisk's writer applies to each
/// tile (`collect()` to R, then back to arrow).
#[derive(Clone, Copy, PartialEq)]
enum Kind {
    /// final type already decided from the schema
    Plain,
    /// uint32 / uint64: int32 when the tile's values fit, else double
    FitUnsigned,
    /// int64: int32 when the tile's values fit, else int64
    FitSigned,
    /// dictionary (R factor): spilled as strings, re-encoded with the levels
    Dict,
}

/// Spill type and kind of a source column. `to_string` columns are cast
/// to string, as the reader's schema does.
fn spill_type(dt: &DataType, to_string: bool) -> Result<(DataType, Kind), String> {
    use DataType::*;
    if to_string {
        return Ok((Utf8, Kind::Plain));
    }
    Ok(match dt {
        Float16 | Float32 | Float64 => (Float64, Kind::Plain),
        Int8 | Int16 | Int32 | UInt8 | UInt16 => (Int32, Kind::Plain),
        UInt32 | UInt64 => (dt.clone(), Kind::FitUnsigned),
        Int64 => (Int64, Kind::FitSigned),
        Utf8 | LargeUtf8 | Utf8View => (Utf8, Kind::Plain),
        Binary | LargeBinary | BinaryView => (Binary, Kind::Plain),
        Boolean => (Boolean, Kind::Plain),
        Date32 => (Date32, Kind::Plain),
        Dictionary(_, v) if matches!(v.as_ref(), Utf8 | LargeUtf8 | Utf8View) => (Utf8, Kind::Dict),
        o => return Err(format!("column type {o} not supported yet")),
    })
}

/// R's integers are 32-bit with -2^31 reserved for NA.
const R_INT_MAX: f64 = 2147483647.0;

/// One tile's column in the type R writes it with.
fn finalize(c: &ArrayRef, kind: Kind, levels: Option<&Levels>) -> Result<ArrayRef, String> {
    let fits = |lo: f64, hi: f64| lo >= -R_INT_MAX && hi <= R_INT_MAX;
    match kind {
        Kind::Plain => Ok(c.clone()),
        Kind::FitUnsigned | Kind::FitSigned => {
            let v = cast(c, &DataType::Float64).map_err(|e| e.to_string())?;
            let v = v.as_any().downcast_ref::<Float64Array>().unwrap();
            let lo = arrow::compute::min(v).unwrap_or(0.0);
            let hi = arrow::compute::max(v).unwrap_or(0.0);
            let to = if fits(lo, hi) {
                DataType::Int32
            } else if kind == Kind::FitSigned {
                DataType::Int64
            } else {
                DataType::Float64
            };
            cast(c, &to).map_err(|e| e.to_string())
        }
        Kind::Dict => levels.ok_or("dictionary levels missing")?.encode(c),
    }
}

/// Factor levels of a dictionary column: the union of the source's
/// dictionaries in scan order (one shared dictionary in the usual case).
#[derive(Default)]
struct Levels {
    values: Vec<String>,
    index: HashMap<String, usize>,
}

impl Levels {
    fn add(&mut self, c: &ArrayRef) -> Result<(), String> {
        let d = c.as_any_dictionary_opt().ok_or("expected a dictionary column")?;
        let v = cast(d.values(), &DataType::Utf8).map_err(|e| e.to_string())?;
        let v = v.as_any().downcast_ref::<StringArray>().unwrap();
        for s in v.iter().flatten() {
            if !self.index.contains_key(s) {
                self.index.insert(s.to_string(), self.values.len());
                self.values.push(s.to_string());
            }
        }
        Ok(())
    }

    /// R's index type: int8 below 127 levels, int16 below 32767, else int32.
    fn encode(&self, c: &ArrayRef) -> Result<ArrayRef, String> {
        let s = c.as_any().downcast_ref::<StringArray>().ok_or("expected strings")?;
        let keys: Vec<Option<i64>> = s.iter()
            .map(|v| v.map(|v| self.index.get(v).map(|k| *k as i64).ok_or(format!("level `{v}` unknown"))).transpose())
            .collect::<Result<_, _>>()?;
        let values: ArrayRef = Arc::new(StringArray::from(self.values.clone()));
        let n = self.values.len();
        let out: ArrayRef = if n < 127 {
            Arc::new(DictionaryArray::<Int8Type>::try_new(
                Int8Array::from(keys.iter().map(|k| k.map(|k| k as i8)).collect::<Vec<_>>()), values)
                .map_err(|e| e.to_string())?)
        } else if n < 32767 {
            Arc::new(DictionaryArray::<Int16Type>::try_new(
                Int16Array::from(keys.iter().map(|k| k.map(|k| k as i16)).collect::<Vec<_>>()), values)
                .map_err(|e| e.to_string())?)
        } else {
            Arc::new(DictionaryArray::<Int32Type>::try_new(
                Int32Array::from(keys.iter().map(|k| k.map(|k| k as i32)).collect::<Vec<_>>()), values)
                .map_err(|e| e.to_string())?)
        };
        Ok(out)
    }
}

pub struct WriteOpts {
    pub id_col: String,
    pub to_string: Vec<String>,
    pub drop: Vec<String>,
    pub root: String,
    pub tile_dirs: Vec<String>,
    pub geo_json: String,
    pub spill_dir: String,
    pub spill_rows: usize,
    /// write the WKB `geom` column and geoparquet metadata (GiottoDisk's
    /// current point format); false keeps only x_index / y_index
    pub geom: bool,
}

/// Routes every row to its tiles, spilling per-tile Arrow IPC parts, then
/// writes each tile as `<root>/<tile_dir>/part-0.parquet`. Returns rows per tile.
pub fn write_tiles(
    spec: &Spec, bounds: Vec<[f64; 4]>, o: &WriteOpts, threads: usize,
) -> Result<Vec<f64>, String> {
    let (fx, fy) = xy_f32(spec)?;
    let tiles = Tiles::new(bounds, fx, fy);
    let nt = tiles.len();
    let src = source_schema(spec)?;
    // meta columns: id first, then the rest in source order, minus x, y, drops
    let mut meta: Vec<String> = vec![o.id_col.clone()];
    for f in src.fields() {
        let n = f.name();
        if n != &o.id_col && n != &spec.x && n != &spec.y && !o.drop.contains(n) {
            meta.push(n.clone());
        }
    }
    let mut meta_fields = Vec::with_capacity(meta.len());
    let mut kinds = Vec::with_capacity(meta.len());
    for m in &meta {
        let f = src.field_with_name(m).map_err(|_| format!("no column `{m}`"))?;
        let (dt, kind) = spill_type(f.data_type(), o.to_string.contains(m))?;
        let name = if m == &o.id_col { "feat_ID".to_string() } else { m.clone() };
        meta_fields.push(Field::new(name, dt, true));
        kinds.push(kind);
    }
    let mut levels: Vec<Option<Levels>> =
        kinds.iter().map(|k| if *k == Kind::Dict { Some(Levels::default()) } else { None }).collect();
    // spill schema: x, y (flipped, double) then meta in output types
    let mut sf = vec![Field::new("__x", DataType::Float64, true), Field::new("__y", DataType::Float64, true)];
    sf.extend(meta_fields.iter().cloned());
    let spill_schema: SchemaRef = Arc::new(Schema::new(sf));
    create_dir_all(&o.spill_dir).map_err(|e| format!("{}: {e}", o.spill_dir))?;

    let mut bufs: Vec<Vec<RecordBatch>> = vec![Vec::new(); nt];
    let mut parts: Vec<usize> = vec![0; nt];
    let mut buffered = 0usize;
    let mut counts = vec![0f64; nt];
    let flush = |bufs: &mut Vec<Vec<RecordBatch>>, parts: &mut Vec<usize>| -> Result<(), String> {
        let jobs: Vec<(usize, Vec<RecordBatch>, usize)> = bufs
            .iter_mut()
            .enumerate()
            .filter(|(_, b)| !b.is_empty())
            .map(|(k, b)| { let p = parts[k]; parts[k] += 1; (k, std::mem::take(b), p) })
            .collect();
        jobs.into_par_iter().try_for_each(|(k, bs, p)| -> Result<(), String> {
            let path = format!("{}/t{k}_{p:06}.arrow", o.spill_dir);
            let f = File::create(&path).map_err(|e| format!("{path}: {e}"))?;
            let mut w = FileWriter::try_new(f, &spill_schema).map_err(|e| e.to_string())?;
            for b in &bs { w.write(b).map_err(|e| e.to_string())?; }
            w.finish().map_err(|e| e.to_string())
        })
    };

    let mut hit = Vec::new();
    for_each_batch(spec, &meta, threads, |b| {
        let (xs, ys) = xy(spec, &b)?;
        // rows per tile, in source order
        let mut idx: HashMap<u32, Vec<u32>> = HashMap::new();
        for r in 0..xs.len() {
            tiles.matches(xs[r], ys[r], &mut hit);
            for k in &hit {
                idx.entry(*k).or_default().push(r as u32);
            }
        }
        // converted columns of the whole batch, then one take per tile
        let mut cols: Vec<ArrayRef> = vec![
            Arc::new(Float64Array::from(xs)), Arc::new(Float64Array::from(ys)),
        ];
        for ((m, f), lv) in meta.iter().zip(&meta_fields).zip(levels.iter_mut()) {
            let c = b.column_by_name(m).ok_or(format!("no column `{m}`"))?;
            if let Some(lv) = lv { lv.add(c)?; }
            cols.push(cast(c, f.data_type()).map_err(|e| format!("cast {m}: {e}"))?);
        }
        let mut keys: Vec<u32> = idx.keys().copied().collect();
        keys.sort_unstable();
        for k in keys {
            let rows = UInt32Array::from(idx.remove(&k).unwrap());
            let taken: Vec<ArrayRef> = cols.iter()
                .map(|c| take(c.as_ref(), &rows, None).map_err(|e| e.to_string()))
                .collect::<Result<_, _>>()?;
            let tb = RecordBatch::try_new(spill_schema.clone(), taken).map_err(|e| e.to_string())?;
            counts[k as usize] += tb.num_rows() as f64;
            buffered += tb.num_rows();
            bufs[k as usize].push(tb);
        }
        if buffered >= o.spill_rows {
            flush(&mut bufs, &mut parts)?;
            buffered = 0;
        }
        Ok(())
    })?;
    flush(&mut bufs, &mut parts)?;

    // each tile: read its parts in order, number, Morton-sort, write
    (0..nt).into_par_iter().filter(|k| parts[*k] > 0).try_for_each(|k| -> Result<(), String> {
        let mut bs = Vec::new();
        for p in 0..parts[k] {
            let path = format!("{}/t{k}_{p:06}.arrow", o.spill_dir);
            let f = File::open(&path).map_err(|e| format!("{path}: {e}"))?;
            for b in FileReader::try_new(f, None).map_err(|e| e.to_string())? {
                bs.push(b.map_err(|e| e.to_string())?);
            }
            std::fs::remove_file(&path).ok();
        }
        let all = concat_batches(&spill_schema, &bs).map_err(|e| e.to_string())?;
        drop(bs);
        write_tile(&all, &meta_fields, &kinds, &levels, &o.root, &o.tile_dirs[k], &o.geo_json, o.geom)
    })?;
    std::fs::remove_dir(&o.spill_dir).ok();
    Ok(counts)
}

#[allow(clippy::too_many_arguments)]
fn write_tile(
    all: &RecordBatch, meta_fields: &[Field], kinds: &[Kind], levels: &[Option<Levels>], root: &str,
    dir: &str, geo: &str, with_geom: bool,
) -> Result<(), String> {
    let x = all.column(0).as_any().downcast_ref::<Float64Array>().unwrap().values();
    let y = all.column(1).as_any().downcast_ref::<Float64Array>().unwrap().values();
    let (ord, _) = order(x, y);
    let (mut fields, mut cols) = base_columns(x, y, &ord, 0.0, with_geom);
    let ordv = UInt32Array::from(ord);
    for (((c, f), kind), lv) in all.columns().iter().skip(2).zip(meta_fields).zip(kinds).zip(levels) {
        let c = finalize(&take(c.as_ref(), &ordv, None).map_err(|e| e.to_string())?, *kind, lv.as_ref())?;
        fields.push(Field::new(f.name(), c.data_type().clone(), true));
        cols.push(c);
    }
    write_parquet(&format!("{root}/{dir}/part-0.parquet"), fields, cols,
        if with_geom { Some(geo.to_string()) } else { None })
}
