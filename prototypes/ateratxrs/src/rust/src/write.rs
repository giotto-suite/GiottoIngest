//! One flush of one store tile -> one parquet part file. A port of
//! GiottoDisk's `.store_write_points_arrow()` for Atera transcripts: Morton
//! row order, WKB point geometry, the same columns and types, geoparquet
//! metadata with the part's bbox.

use std::collections::HashMap;
use std::fs::{create_dir_all, File};
use std::sync::Arc;

use arrow::array::{
    ArrayRef, BinaryBuilder, Float32Array, Float64Array, Int32Array, StringArray, StringBuilder,
    UInt8Array,
};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;

use crate::cell::{Buf, Lookups};

const UNASSIGNED: u32 = 4_294_967_295;

/// GiottoDisk's `.morton_encode()`, including R's signed 32-bit wrap; `None`
/// where R gives NA (a zero-width extent).
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
    let xi = q(x, e[0], e[1])?;
    let yi = q(y, e[2], e[3])?;
    Some((spread(xi) | (spread(yi) << 1)) as i32)
}

/// `.encode_xenium_id()`: 8 hex digits mapped onto a-p, then `-suffix`.
fn cell_id(prefix: u32, suffix: u32, s: &mut String) {
    s.clear();
    if prefix == UNASSIGNED {
        s.push_str("UNASSIGNED");
        return;
    }
    for d in format!("{:04x}{:04x}", prefix / 65536, prefix % 65536).bytes() {
        let k = if d.is_ascii_digit() { d - b'0' } else { d - b'a' + 10 };
        s.push((b'a' + k) as char);
    }
    s.push('-');
    if suffix > i32::MAX as u32 {
        s.push_str("NA");
    } else {
        s.push_str(&suffix.to_string());
    }
}

fn num4(v: f64) -> String {
    let s = format!("{v:.4}");
    let s = if s.contains('.') { s.trim_end_matches('0').trim_end_matches('.').to_string() } else { s };
    if s == "-0" { "0".into() } else { s }
}

pub fn write_points(
    dir: &str, n_part: usize, row_offset: f64, b: &Buf, ord: u64, lk: &Lookups,
) -> Result<(), String> {
    let n = b.len();
    let fold = |v: &[f64], f: fn(f64, f64) -> f64, init: f64| v.iter().fold(init, |a, c| f(a, *c));
    let gext = [
        fold(&b.x, f64::min, f64::INFINITY), fold(&b.x, f64::max, f64::NEG_INFINITY),
        fold(&b.y, f64::min, f64::INFINITY), fold(&b.y, f64::max, f64::NEG_INFINITY),
    ];
    let mut o: Vec<usize> = (0..n).collect();
    let keys: Vec<Option<i32>> = (0..n).map(|i| morton(b.x[i], b.y[i], gext)).collect();
    // stable, NA last: R's order()
    o.sort_by_key(|&i| (keys[i].is_none(), keys[i].unwrap_or(0)));

    let row_index = Int32Array::from_iter_values(o.iter().map(|&i| (row_offset + (i + 1) as f64) as i32));
    let x_index = Float64Array::from_iter_values(o.iter().map(|&i| b.x[i]));
    let y_index = Float64Array::from_iter_values(o.iter().map(|&i| b.y[i]));
    let mut geom = BinaryBuilder::with_capacity(n, n * 21);
    let mut wkb = [0u8; 21];
    wkb[0] = 1;
    wkb[1..5].copy_from_slice(&1u32.to_le_bytes());
    for &i in &o {
        wkb[5..13].copy_from_slice(&b.x[i].to_le_bytes());
        wkb[13..21].copy_from_slice(&b.y[i].to_le_bytes());
        geom.append_value(wkb);
    }
    let mut feat = StringBuilder::with_capacity(n, n * 8);
    for &i in &o {
        let g = usize::try_from(b.cw[i]).ok()
            .and_then(|c| lk.cw_gene.get(c).copied().flatten())
            .and_then(|k| usize::try_from(k).ok())
            .and_then(|k| lk.gene_names.get(k));
        match g { Some(s) => feat.append_value(s), None => feat.append_null() }
    }
    let mut tid = StringBuilder::with_capacity(n, n * 12);
    for &i in &o {
        tid.append_value(((ord << 32) + b.row[i]).to_string());
    }
    let mut cid = StringBuilder::with_capacity(n, n * 11);
    let mut s = String::with_capacity(16);
    for &i in &o {
        cell_id(b.cid0[i], b.cid1[i], &mut s);
        cid.append_value(&s);
    }
    let ovl = UInt8Array::from_iter_values(o.iter().map(|&i| b.ovl[i]));
    let z = Float32Array::from_iter_values(o.iter().map(|&i| b.z[i]));
    let qv = Float32Array::from_iter_values(o.iter().map(|&i| b.qv[i]));
    let fov = StringArray::new_null(n);
    let nd = Float32Array::new_null(n);
    let cw = Int32Array::from_iter_values(o.iter().map(|&i| b.cw[i]));

    let geo = format!(
        "{{\"version\":\"1.0.0\",\"primary_column\":\"geom\",\"columns\":{{\"geom\":{{\"encoding\":\"WKB\",\"geometry_types\":[\"Point\"],\"crs\":null,\"edges\":\"planar\",\"bbox\":[{},{},{},{}]}}}}}}",
        num4(gext[0]), num4(gext[2]), num4(gext[1]), num4(gext[3])
    );
    let f = |nm: &str, dt: DataType| Field::new(nm, dt, true);
    let schema = Arc::new(Schema::new_with_metadata(
        vec![
            f("row_index", DataType::Int32), f("x_index", DataType::Float64),
            f("y_index", DataType::Float64), f("geom", DataType::Binary),
            f("feat_ID", DataType::Utf8), f("transcript_id", DataType::Utf8),
            f("cell_id", DataType::Utf8), f("overlaps_nucleus", DataType::UInt8),
            f("z_location", DataType::Float32), f("qv", DataType::Float32),
            f("fov_name", DataType::Utf8), f("nucleus_distance", DataType::Float32),
            f("codeword_index", DataType::Int32),
        ],
        HashMap::from([("geo".to_string(), geo)]),
    ));
    let cols: Vec<ArrayRef> = vec![
        Arc::new(row_index), Arc::new(x_index), Arc::new(y_index), Arc::new(geom.finish()),
        Arc::new(feat.finish()), Arc::new(tid.finish()), Arc::new(cid.finish()), Arc::new(ovl),
        Arc::new(z), Arc::new(qv), Arc::new(fov), Arc::new(nd), Arc::new(cw),
    ];
    let batch = RecordBatch::try_new(schema.clone(), cols).map_err(|e| e.to_string())?;

    create_dir_all(dir).map_err(|e| format!("{dir}: {e}"))?;
    let path = format!("{dir}/part-{n_part:05}-0.parquet");
    let file = File::create(&path).map_err(|e| format!("{path}: {e}"))?;
    let props = WriterProperties::builder()
        .set_compression(Compression::SNAPPY)
        .set_max_row_group_size(1 << 20)
        .build();
    let mut w = ArrowWriter::try_new(file, schema, Some(props)).map_err(|e| e.to_string())?;
    w.write(&batch).map_err(|e| e.to_string())?;
    w.close().map_err(|e| e.to_string())?;
    Ok(())
}
