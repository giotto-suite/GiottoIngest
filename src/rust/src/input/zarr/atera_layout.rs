//! One flush of one store tile -> one parquet part file: GiottoDisk's
//! `.store_write_points_arrow()` for Atera transcripts. The Xenium-layout
//! metadata columns are built here, in Morton order; the row order, leading
//! columns and writer are the shared core's.

use std::sync::Arc;

use arrow::array::{ArrayRef, Float32Array, Int32Array, StringArray, StringBuilder, UInt8Array};
use arrow::datatypes::{DataType, Field};

use super::atera::{Buf, Lookups};
use crate::core::geo::points_geo;
use crate::core::morton::order;
use crate::core::tile_write::{base_columns, write_parquet};

const UNASSIGNED: u32 = 4_294_967_295;

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

pub fn write_points(
    dir: &str, n_part: usize, row_offset: f64, b: &Buf, ord: u64, lk: &Lookups,
) -> Result<(), String> {
    let n = b.len();
    let (ord32, gext) = order(&b.x, &b.y);
    let (mut fields, mut cols) = base_columns(&b.x, &b.y, &ord32, row_offset, true);
    let o: Vec<usize> = ord32.iter().map(|&i| i as usize).collect();

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
    let f = |nm: &str, dt: DataType| Field::new(nm, dt, true);
    fields.extend([
        f("feat_ID", DataType::Utf8), f("transcript_id", DataType::Utf8),
        f("cell_id", DataType::Utf8), f("overlaps_nucleus", DataType::UInt8),
        f("z_location", DataType::Float32), f("qv", DataType::Float32),
        f("fov_name", DataType::Utf8), f("nucleus_distance", DataType::Float32),
        f("codeword_index", DataType::Int32),
    ]);
    let meta: [ArrayRef; 9] = [
        Arc::new(feat.finish()), Arc::new(tid.finish()), Arc::new(cid.finish()),
        Arc::new(UInt8Array::from_iter_values(o.iter().map(|&i| b.ovl[i]))),
        Arc::new(Float32Array::from_iter_values(o.iter().map(|&i| b.z[i]))),
        Arc::new(Float32Array::from_iter_values(o.iter().map(|&i| b.qv[i]))),
        Arc::new(StringArray::new_null(n)), Arc::new(Float32Array::new_null(n)),
        Arc::new(Int32Array::from_iter_values(o.iter().map(|&i| b.cw[i]))),
    ];
    cols.extend(meta);

    // the part's own bbox; R writes 0 where the flipped extent is -0
    let z = |v: f64| if v == 0.0 { 0.0 } else { v };
    let geo = points_geo([z(gext[0]), z(gext[1]), z(gext[2]), z(gext[3])]);
    write_parquet(&format!("{dir}/part-{n_part:05}-0.parquet"), fields, cols, Some(geo))
}
