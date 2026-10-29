//! Point tile store passes: extent, per-depth tile counts, and the routed
//! write. The write reproduces GiottoDisk's per-tile path
//! (`storeWrite(parquetGeomStore, data.frame)` -> `.store_write_points`):
//! rows in source order numbered 1..n, Morton order, WKB points, and the
//! column types an R round trip produces.

use std::collections::HashMap;
use std::fs::{create_dir_all, File};
use std::sync::Arc;

use arrow::array::{
    Array, ArrayRef, BinaryBuilder, Float64Array, Int32Array, RecordBatch, UInt32Array,
};
use arrow::compute::{cast, concat_batches, take};
use arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use arrow::ipc::reader::FileReader;
use arrow::ipc::writer::FileWriter;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;
use rayon::prelude::*;

use crate::source::{for_each_batch, reduce_units, source_schema, xy, xy_f32, Spec};
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

/// Output type of a source column after R's collect / as.data.table / arrow
/// round trip, as GiottoDisk writes it.
fn r_type(dt: &DataType, to_string: bool) -> Result<DataType, String> {
    if to_string {
        return Ok(DataType::Utf8);
    }
    Ok(match dt {
        DataType::Float16 | DataType::Float32 | DataType::Float64 => DataType::Float64,
        DataType::Int8 | DataType::Int16 | DataType::Int32 | DataType::UInt8 | DataType::UInt16 => {
            DataType::Int32
        }
        DataType::UInt32 => DataType::Float64,
        DataType::Utf8 | DataType::LargeUtf8 | DataType::Utf8View => DataType::Utf8,
        DataType::Boolean => DataType::Boolean,
        o => return Err(format!("column type {o} not supported yet")),
    })
}

/// GiottoDisk's `.morton_encode()`, with R's signed 32-bit wrap; `None` is NA.
fn morton(x: f64, y: f64, e: [f64; 4]) -> Option<i32> {
    let q = |v: f64, lo: f64, hi: f64| -> Option<u32> {
        let t = (v - lo) / (hi - lo) * 65535.0;
        if t.is_finite() { Some(t.trunc() as i64 as u32) } else { None }
    };
    let spread = |v: u32| -> u32 {
        let mut v = v & 0x0000_FFFF;
        v = (v | (v << 8)) & 0x00FF_00FF;
        v = (v | (v << 4)) & 0x0F0F_0F0F;
        v = (v | (v << 2)) & 0x3333_3333;
        (v | (v << 1)) & 0x5555_5555
    };
    Some((spread(q(x, e[0], e[1])?) | (spread(q(y, e[2], e[3])?) << 1)) as i32)
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
    for m in &meta {
        let f = src.field_with_name(m).map_err(|_| format!("no column `{m}`"))?;
        let dt = r_type(f.data_type(), o.to_string.contains(m))?;
        let name = if m == &o.id_col { "feat_ID".to_string() } else { m.clone() };
        meta_fields.push(Field::new(name, dt, true));
    }
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
        for (m, f) in meta.iter().zip(&meta_fields) {
            let c = b.column_by_name(m).ok_or(format!("no column `{m}`"))?;
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
        write_tile(&all, &meta_fields, &o.root, &o.tile_dirs[k], &o.geo_json, o.geom)
    })?;
    std::fs::remove_dir(&o.spill_dir).ok();
    Ok(counts)
}

fn write_tile(
    all: &RecordBatch, meta_fields: &[Field], root: &str, dir: &str, geo: &str, with_geom: bool,
) -> Result<(), String> {
    let n = all.num_rows();
    let x = all.column(0).as_any().downcast_ref::<Float64Array>().unwrap();
    let y = all.column(1).as_any().downcast_ref::<Float64Array>().unwrap();
    let fold = |a: &Float64Array, f: fn(f64, f64) -> f64, i: f64| a.values().iter().fold(i, |p, v| f(p, *v));
    let gext = [
        fold(x, f64::min, f64::INFINITY), fold(x, f64::max, f64::NEG_INFINITY),
        fold(y, f64::min, f64::INFINITY), fold(y, f64::max, f64::NEG_INFINITY),
    ];
    let keys: Vec<Option<i32>> = (0..n).map(|i| morton(x.value(i), y.value(i), gext)).collect();
    let mut ord: Vec<u32> = (0..n as u32).collect();
    ord.sort_by_key(|&i| (keys[i as usize].is_none(), keys[i as usize].unwrap_or(0)));
    let ordv = UInt32Array::from(ord.clone());

    let row_index = Int32Array::from_iter_values(ord.iter().map(|&i| i as i32 + 1));
    let xo = Float64Array::from_iter_values(ord.iter().map(|&i| x.value(i as usize)));
    let yo = Float64Array::from_iter_values(ord.iter().map(|&i| y.value(i as usize)));
    let mut fields = vec![
        Field::new("row_index", DataType::Int32, true), Field::new("x_index", DataType::Float64, true),
        Field::new("y_index", DataType::Float64, true),
    ];
    let mut cols: Vec<ArrayRef> = vec![Arc::new(row_index), Arc::new(xo), Arc::new(yo)];
    if with_geom {
        let mut geom = BinaryBuilder::with_capacity(n, n * 21);
        let mut wkb = [0u8; 21];
        wkb[0] = 1;
        wkb[1..5].copy_from_slice(&1u32.to_le_bytes());
        for &i in &ord {
            wkb[5..13].copy_from_slice(&x.value(i as usize).to_le_bytes());
            wkb[13..21].copy_from_slice(&y.value(i as usize).to_le_bytes());
            geom.append_value(wkb);
        }
        fields.push(Field::new("geom", DataType::Binary, true));
        cols.push(Arc::new(geom.finish()));
    }
    fields.extend(meta_fields.iter().cloned());
    for c in all.columns().iter().skip(2) {
        cols.push(take(c.as_ref(), &ordv, None).map_err(|e| e.to_string())?);
    }
    let md: HashMap<String, String> = if with_geom {
        HashMap::from([("geo".to_string(), geo.to_string())])
    } else {
        HashMap::new()
    };
    let schema = Arc::new(Schema::new_with_metadata(fields, md));
    let batch = RecordBatch::try_new(schema.clone(), cols).map_err(|e| e.to_string())?;

    let d = format!("{root}/{dir}");
    create_dir_all(&d).map_err(|e| format!("{d}: {e}"))?;
    let path = format!("{d}/part-0.parquet");
    let f = File::create(&path).map_err(|e| format!("{path}: {e}"))?;
    let props = WriterProperties::builder()
        .set_compression(Compression::SNAPPY)
        .set_max_row_group_size(1 << 20)
        .build();
    let mut w = ArrowWriter::try_new(f, schema, Some(props)).map_err(|e| e.to_string())?;
    w.write(&batch).map_err(|e| e.to_string())?;
    w.close().map_err(|e| e.to_string())?;
    Ok(())
}
