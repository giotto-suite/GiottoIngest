//! Polygon tile routing, as GiottoDisk's polygon tile write selects rows:
//! each polygon (group) is assigned to the tile holding its envelope centroid
//! (`.dplyr_xy_envelopes()` + `.dplyr_crop()` with the contiguous `getTile()`
//! edge rules over the written tiles), and all its vertex rows go there in
//! `row_index` order. Building the polygons stays in R.

use std::collections::HashMap;
use std::fs::{create_dir_all, File};

use arrow::array::{Array, ArrayRef, StringArray, UInt32Array};
use arrow::compute::{cast, concat_batches, sort_to_indices, take};
use arrow::datatypes::DataType;
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;
use rayon::prelude::*;

use crate::source::{for_each_batch, reduce_units, source_schema, Spec};
use crate::tiles::Tiles;

/// Envelope per polygon id; `None` is the null id, which R groups as one
/// polygon (`group_by()` keeps NA as a group).
type Env = HashMap<Option<String>, [f64; 4]>;

fn key(k: &StringArray, i: usize) -> Option<String> {
    if k.is_null(i) { None } else { Some(k.value(i).to_string()) }
}

fn keys(b: &RecordBatch, group: &str) -> Result<StringArray, String> {
    let c = b.column_by_name(group).ok_or(format!("no column `{group}`"))?;
    let s = cast(c, &DataType::Utf8).map_err(|e| e.to_string())?;
    Ok(s.as_any().downcast_ref::<StringArray>().unwrap().clone())
}

fn coords(spec: &Spec, b: &RecordBatch) -> Result<(Vec<f64>, Vec<f64>), String> {
    crate::source::xy(spec, b)
}

pub struct Routed {
    pub counts: Vec<f64>,
    pub groups: f64,
    pub unassigned: f64,
}

pub fn route_polygons(
    spec: &Spec, bounds: Vec<[f64; 4]>, group: &str, out_dir: &str, tile_files: &[String],
) -> Result<Routed, String> {
    let sch = source_schema(spec)?;
    let x32 = sch.field_with_name(&spec.x).map_err(|e| e.to_string())?.data_type() == &DataType::Float32;
    let y32 = sch.field_with_name(&spec.y).map_err(|e| e.to_string())?.data_type() == &DataType::Float32;

    // 1. envelope per group (min / max of the, possibly flipped, vertices)
    let cols = vec![group.to_string()];
    let env: Env = reduce_units(spec, &cols, Env::new, |acc, b| {
        let k = keys(&b, group)?;
        let (xs, ys) = coords(spec, &b)?;
        for i in 0..xs.len() {
            let e = acc.entry(key(&k, i))
                .or_insert([f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY]);
            if !xs[i].is_nan() { e[0] = e[0].min(xs[i]); e[1] = e[1].max(xs[i]); }
            if !ys[i].is_nan() { e[2] = e[2].min(ys[i]); e[3] = e[3].max(ys[i]); }
        }
        Ok(())
    }, |mut a, b| {
        for (k, v) in b {
            let e = a.entry(k).or_insert(v);
            e[0] = e[0].min(v[0]); e[1] = e[1].max(v[1]); e[2] = e[2].min(v[2]); e[3] = e[3].max(v[3]);
        }
        a
    })?;

    // 2. centroid -> tile. arrow adds in the column type, then divides in double
    let tiles = Tiles::new(bounds, false, false);
    let nt = tiles.len();
    let mid = |lo: f64, hi: f64, f32col: bool| {
        if f32col { ((lo as f32) + (hi as f32)) as f64 / 2.0 } else { (lo + hi) / 2.0 }
    };
    let mut hit = Vec::new();
    let mut tile_of: HashMap<String, u32> = HashMap::with_capacity(env.len());
    let mut null_tile: Option<u32> = None;
    let mut unassigned = 0f64;
    for (k, e) in &env {
        tiles.matches(mid(e[0], e[1], x32), mid(e[2], e[3], y32), &mut hit);
        if hit.is_empty() { unassigned += 1.0; }
        // a centroid on no tile edge meets one tile at most
        if let Some(t) = hit.first() {
            match k {
                Some(k) => { tile_of.insert(k.clone(), *t); }
                None => null_tile = Some(*t),
            }
        }
    }
    let groups = env.len() as f64;
    drop(env);

    // 3. route every vertex row of its polygon's tile, all source columns,
    //    y as the query returns it (flipped)
    let all_cols: Vec<String> = sch.fields().iter().map(|f| f.name().clone()).collect();
    let mut bufs: Vec<Vec<RecordBatch>> = vec![Vec::new(); nt];
    for_each_batch(spec, &all_cols, rayon::current_num_threads(), |b| {
        let k = keys(&b, group)?;
        let mut idx: HashMap<u32, Vec<u32>> = HashMap::new();
        for i in 0..b.num_rows() {
            let t = if k.is_null(i) { null_tile.as_ref() } else { tile_of.get(k.value(i)) };
            if let Some(t) = t {
                idx.entry(*t).or_default().push(i as u32);
            }
        }
        let b = if spec.flip_y { flip(&b, &spec.y)? } else { b };
        for (t, rows) in idx {
            let r = UInt32Array::from(rows);
            let cols: Vec<ArrayRef> = b.columns().iter()
                .map(|c| take(c.as_ref(), &r, None).map_err(|e| e.to_string()))
                .collect::<Result<_, _>>()?;
            bufs[t as usize].push(RecordBatch::try_new(b.schema(), cols).map_err(|e| e.to_string())?);
        }
        Ok(())
    })?;

    // 4. per tile: row_index order, one parquet file
    create_dir_all(out_dir).map_err(|e| format!("{out_dir}: {e}"))?;
    let counts: Vec<f64> = bufs.into_par_iter().enumerate().map(|(t, bs)| -> Result<f64, String> {
        if bs.is_empty() { return Ok(0.0); }
        let all = concat_batches(&bs[0].schema(), &bs).map_err(|e| e.to_string())?;
        let all = match all.column_by_name("row_index") {
            Some(ri) => {
                let ord = sort_to_indices(ri.as_ref(), None, None).map_err(|e| e.to_string())?;
                let cols: Vec<ArrayRef> = all.columns().iter()
                    .map(|c| take(c.as_ref(), &ord, None).map_err(|e| e.to_string()))
                    .collect::<Result<_, _>>()?;
                RecordBatch::try_new(all.schema(), cols).map_err(|e| e.to_string())?
            }
            None => all,
        };
        let path = format!("{out_dir}/{}", tile_files[t]);
        let f = File::create(&path).map_err(|e| format!("{path}: {e}"))?;
        let props = WriterProperties::builder().set_compression(Compression::SNAPPY).build();
        let mut w = ArrowWriter::try_new(f, all.schema(), Some(props)).map_err(|e| e.to_string())?;
        w.write(&all).map_err(|e| e.to_string())?;
        w.close().map_err(|e| e.to_string())?;
        Ok(all.num_rows() as f64)
    }).collect::<Result<_, _>>()?;
    Ok(Routed { counts, groups, unassigned })
}

/// Negate the y column in its own type, as `mutate(y = -y)` does.
fn flip(b: &RecordBatch, y: &str) -> Result<RecordBatch, String> {
    let i = b.schema().index_of(y).map_err(|e| e.to_string())?;
    let c = b.column(i);
    let neg = arrow::compute::kernels::numeric::neg(c.as_ref()).map_err(|e| e.to_string())?;
    let mut cols = b.columns().to_vec();
    cols[i] = neg;
    RecordBatch::try_new(b.schema(), cols).map_err(|e| e.to_string())
}

