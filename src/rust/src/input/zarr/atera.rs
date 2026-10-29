//! One Atera grid cell -> its tile partitions. A port of GiottoDisk's
//! `.atera_tx_cell_worker()`: same sample, same quantile cuts, same row
//! routing, buffering and flush order, so the written rows and their
//! `row_index` match the R decoder.

use super::archive::{f16_le, f32_le, u32_le, u8_le, Archive};
use super::atera_layout::write_points;

pub struct CellPlan {
    pub tile: String,
    pub ord: u64,
    pub s: usize,
    pub tile_index: Vec<i32>,
}

pub struct Opts {
    pub root: String,
    pub qv_threshold: Option<f64>,
    pub flip: bool,
    pub window_rows: usize,
    pub buffer_rows: usize,
}

pub struct Lookups {
    pub gene_names: Vec<String>,
    pub cw_gene: Vec<Option<i32>>,
}

pub struct CellResult {
    pub tile_index: Vec<i32>,
    pub counts: Vec<f64>,
    pub bounds: Vec<[f64; 4]>,
}

/// Buffered rows of one store tile, in arrival order.
#[derive(Default)]
pub struct Buf {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub z: Vec<f32>,
    pub qv: Vec<f32>,
    pub cw: Vec<i32>,
    pub cid0: Vec<u32>,
    pub cid1: Vec<u32>,
    pub ovl: Vec<u8>,
    pub row: Vec<u64>,
}

impl Buf {
    pub fn len(&self) -> usize { self.x.len() }
    fn clear(&mut self) { *self = Buf::default(); }
}

/// R's `quantile(v, seq_len(s - 1) / s, type = 1)`, then `unique()`.
fn quantile_cuts(v: &[f64], s: usize) -> Vec<f64> {
    if s < 2 || v.is_empty() {
        return Vec::new();
    }
    let n = v.len();
    let mut sorted = v.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let fuzz = 4.0 * f64::EPSILON;
    // R pads the sorted values as c(x1, x1, x, xn, xn); index is 1-based
    let pad = |i: usize| -> f64 {
        if i <= 2 { sorted[0] } else if i >= n + 3 { sorted[n - 1] } else { sorted[i - 3] }
    };
    let mut out = Vec::with_capacity(s - 1);
    for k in 1..s {
        let p = k as f64 / s as f64;
        let nppm = n as f64 * p;
        let j = (nppm + fuzz).floor();
        let h = nppm > j;
        let j = j as usize;
        out.push(if h { pad(j + 3) } else { pad(j + 2) });
    }
    out.dedup();
    out
}

/// `findInterval(v, cuts)`: number of cuts <= v.
fn find_interval(cuts: &[f64], v: f64) -> usize {
    cuts.partition_point(|c| *c <= v)
}

struct Split {
    xc: Vec<f64>,
    yc: Vec<Vec<f64>>,
    s: usize,
}

fn quantile_split(x: &[f64], y: &[f64], s: usize) -> Split {
    let xc = quantile_cuts(x, s);
    let yc = (1..=s)
        .map(|b| {
            let yb: Vec<f64> = x.iter().zip(y)
                .filter(|(xv, _)| find_interval(&xc, **xv) + 1 == b)
                .map(|(_, yv)| *yv)
                .collect();
            quantile_cuts(&yb, s)
        })
        .collect();
    Split { xc, yc, s }
}

/// 1-based slot, band-major.
fn quantile_slot(sp: &Split, x: f64, y: f64) -> usize {
    let band = find_interval(&sp.xc, x) + 1;
    let cell = find_interval(&sp.yc[band - 1], y) + 1;
    (band - 1) * sp.s + cell
}

fn quantile_bounds(sp: &Split, e: [f64; 4]) -> Vec<[f64; 4]> {
    let s = sp.s;
    let mut out = vec![[f64::NAN; 4]; s * s];
    let mut xb = vec![e[0]];
    xb.extend_from_slice(&sp.xc);
    xb.push(e[1]);
    for b in 1..=s {
        if b > xb.len() - 1 {
            continue;
        }
        let mut yb = vec![e[2]];
        yb.extend_from_slice(&sp.yc[b - 1]);
        yb.push(e[3]);
        for k in 1..yb.len() {
            out[(b - 1) * s + k - 1] = [xb[b - 1], xb[b], yb[k - 1], yb[k]];
        }
    }
    out
}

/// R's `floor(seq(1, to, length.out = len))`, as 0-based starts.
fn strided_starts(to: f64, len: usize) -> Vec<usize> {
    let by = (to - 1.0) / (len as f64 - 1.0);
    (0..len)
        .map(|k| {
            let v = if k == 0 { 1.0 } else if k == len - 1 { to } else { 1.0 + k as f64 * by };
            v.floor() as usize - 1
        })
        .collect()
}

pub fn run_cell(zarr: &str, cp: &CellPlan, lk: &Lookups, o: &Opts) -> Result<CellResult, String> {
    let mut ar = Archive::open(zarr)?;
    let pref = |a: &str| format!("grid/{}/{}", cp.tile, a);
    let m_loc = ar.meta(&pref("location"))?;
    let m_cw = ar.meta(&pref("codeword_identity"))?;
    let m_qv = ar.meta(&pref("quality_score"))?;
    let m_cid = ar.meta(&pref("cell_id"))?;
    let m_ovl = ar.meta(&pref("overlaps_nucleus"))?;
    for (m, want) in [(&m_loc, "<f4"), (&m_cw, "<u4"), (&m_qv, "<f2"), (&m_cid, "<u4"), (&m_ovl, "|u1")] {
        if m.dtype != want {
            return Err(format!("tile {}: dtype {} where {} expected", cp.tile, m.dtype, want));
        }
    }
    let n = m_loc.shape[0];

    // sample: the whole cell when small, else 16 evenly spaced runs
    let sample_rows = 1_048_576usize;
    let (sx, sy) = if n <= sample_rows {
        let loc = ar.read_cols(&pref("location"), &m_loc, 0, n, 4, f32_le)?;
        (loc[0].iter().map(|v| *v as f64).collect::<Vec<_>>(),
         loc[1].iter().map(|v| *v as f64).collect::<Vec<_>>())
    } else {
        let len = sample_rows / 16;
        let (mut sx, mut sy) = (Vec::new(), Vec::new());
        for st in strided_starts((n - len + 1) as f64, 16) {
            let loc = ar.read_cols(&pref("location"), &m_loc, st, st + len, 4, f32_le)?;
            sx.extend(loc[0].iter().map(|v| *v as f64));
            sy.extend(loc[1].iter().map(|v| *v as f64));
        }
        (sx, sy)
    };
    let sp = quantile_split(&sx, &sy, cp.s);
    let fold = |v: &[f64], f: fn(f64, f64) -> f64, init: f64| v.iter().fold(init, |a, b| f(a, *b));
    let mut bounds = quantile_bounds(&sp, [
        fold(&sx, f64::min, f64::INFINITY), fold(&sx, f64::max, f64::NEG_INFINITY),
        fold(&sy, f64::min, f64::INFINITY), fold(&sy, f64::max, f64::NEG_INFINITY),
    ]);

    let nslot = cp.tile_index.len();
    let mut counts = vec![0f64; nslot];
    let mut ext = vec![[f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY]; nslot];
    let mut buf: Vec<Buf> = (0..nslot).map(|_| Buf::default()).collect();
    let mut buffered = 0usize;
    let mut n_part = 0usize;

    let flush = |buf: &mut Vec<Buf>, counts: &mut Vec<f64>, ext: &mut Vec<[f64; 4]>,
                     n_part: &mut usize| -> Result<(), String> {
        for i in 0..nslot {
            if buf[i].len() == 0 {
                continue;
            }
            *n_part += 1;
            let b = &mut buf[i];
            for (xv, yv) in b.x.iter().zip(&b.y) {
                let e = &mut ext[i];
                e[0] = e[0].min(*xv); e[1] = e[1].max(*xv);
                e[2] = e[2].min(*yv); e[3] = e[3].max(*yv);
            }
            if o.flip {
                for yv in b.y.iter_mut() { *yv = -*yv; }
            }
            let dir = format!("{}/tile_index={:03}", o.root, cp.tile_index[i]);
            write_points(&dir, *n_part, counts[i], b, cp.ord, lk)?;
            counts[i] += b.len() as f64;
            b.clear();
        }
        Ok(())
    };

    let mut r0 = 0usize;
    while r0 < n {
        let r1 = (r0 + o.window_rows).min(n);
        let loc = ar.read_cols(&pref("location"), &m_loc, r0, r1, 4, f32_le)?;
        let cw = ar.read_cols(&pref("codeword_identity"), &m_cw, r0, r1, 4, u32_le)?;
        let qv = ar.read_cols(&pref("quality_score"), &m_qv, r0, r1, 2, f16_le)?;
        let cid = ar.read_cols(&pref("cell_id"), &m_cid, r0, r1, 4, u32_le)?;
        let ovl = ar.read_cols(&pref("overlaps_nucleus"), &m_ovl, r0, r1, 1, u8_le)?;
        let w = r1 - r0;

        // slot per row; rows below the quality threshold are dropped
        let mut keyed: Vec<(usize, usize)> = Vec::with_capacity(w);
        for r in 0..w {
            let q = qv[0][r] as f64;
            if let Some(t) = o.qv_threshold {
                if q < t { continue; }
            }
            let key = quantile_slot(&sp, loc[0][r] as f64, loc[1][r] as f64);
            keyed.push((key, r));
        }
        if keyed.is_empty() {
            r0 = r1;
            continue;
        }
        keyed.sort_by_key(|k| k.0); // stable, as R's radix order()
        for (key, r) in &keyed {
            let b = &mut buf[key - 1];
            b.x.push(loc[0][*r] as f64);
            b.y.push(loc[1][*r] as f64);
            b.z.push(loc[2][*r]);
            b.qv.push(qv[0][*r]);
            b.cw.push(cw[0][*r] as i32);
            b.cid0.push(cid[0][*r]);
            b.cid1.push(cid[1][*r]);
            b.ovl.push(ovl[0][*r]);
            b.row.push((r0 + r) as u64);
        }
        buffered += keyed.len();
        if buffered >= o.buffer_rows {
            flush(&mut buf, &mut counts, &mut ext, &mut n_part)?;
            buffered = 0;
        }
        r0 = r1;
    }
    flush(&mut buf, &mut counts, &mut ext, &mut n_part)?;

    for i in 0..nslot {
        if counts[i] > 0.0 {
            // R's pmin/pmax: NA stays NA
            let na = |f: fn(f64, f64) -> f64, a: f64, c: f64| if a.is_nan() { a } else { f(a, c) };
            let b = &mut bounds[i];
            b[0] = na(f64::min, b[0], ext[i][0]); b[1] = na(f64::max, b[1], ext[i][1]);
            b[2] = na(f64::min, b[2], ext[i][2]); b[3] = na(f64::max, b[3], ext[i][3]);
        }
        if o.flip {
            let b = &mut bounds[i];
            let (lo, hi) = (b[2], b[3]);
            b[2] = -hi;
            b[3] = -lo;
        }
    }
    Ok(CellResult { tile_index: cp.tile_index.clone(), counts, bounds })
}
