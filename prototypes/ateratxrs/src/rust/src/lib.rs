use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

use extendr_api::prelude::*;
use extendr_api::Result;

mod cell;
mod write;
mod zarr;

use crate::cell::{run_cell, CellPlan, CellResult, Lookups, Opts};

/// Stream Atera grid cells into tile partitions on `threads` threads.
///
/// Internal; called by [atera_tx_cells()], which validates the arguments.
/// Each thread takes the next grid cell in the given order, so pass them
/// largest first.
/// @keywords internal
#[extendr]
fn atera_tx_cells_rs(
    zarr_path: &str,
    tiles: Strings,
    ords: Robj,
    s: Robj,
    tile_index: List,
    gene_names: Strings,
    cw_gene: Robj,
    root: &str,
    qv_threshold: f64,
    flip: bool,
    window_rows: f64,
    buffer_rows: f64,
    threads: i32,
) -> Result<List> {
    let ints = |r: &Robj, what: &str| -> Result<Vec<i32>> {
        r.as_integer_slice().map(|v| v.to_vec())
            .ok_or_else(|| Error::Other(format!("`{what}` must be an integer vector")))
    };
    let ords = ints(&ords, "ords")?;
    let s = ints(&s, "s")?;
    let mut plans = Vec::with_capacity(tiles.len());
    for (k, t) in tiles.iter().enumerate() {
        let ti = tile_index.elt(k)?;
        plans.push(CellPlan {
            tile: t.to_string(),
            ord: ords[k] as u64,
            s: s[k] as usize,
            tile_index: ints(&ti, "tile_index")?,
        });
    }
    let lk = Lookups {
        gene_names: gene_names.iter().map(|g| g.to_string()).collect(),
        cw_gene: ints(&cw_gene, "cw_gene")?
            .into_iter()
            .map(|v| if v == i32::MIN { None } else { Some(v) })
            .collect(),
    };
    let opts = Opts {
        root: root.to_string(),
        qv_threshold: if qv_threshold.is_nan() { None } else { Some(qv_threshold) },
        flip,
        window_rows: window_rows as usize,
        buffer_rows: buffer_rows as usize,
    };

    let next = AtomicUsize::new(0);
    let results: Mutex<Vec<Option<std::result::Result<CellResult, String>>>> =
        Mutex::new((0..plans.len()).map(|_| None).collect());
    std::thread::scope(|sc| {
        for _ in 0..threads.max(1) {
            sc.spawn(|| loop {
                let k = next.fetch_add(1, Ordering::SeqCst);
                if k >= plans.len() {
                    break;
                }
                let r = run_cell(zarr_path, &plans[k], &lk, &opts);
                results.lock().unwrap()[k] = Some(r);
            });
        }
    });

    let mut out: Vec<Robj> = Vec::with_capacity(plans.len());
    for (k, r) in results.into_inner().unwrap().into_iter().enumerate() {
        let r = r.ok_or_else(|| Error::Other(format!("cell {} not run", plans[k].tile)))?
            .map_err(|e| Error::Other(format!("[ateratxrs] cell {}: {e}", plans[k].tile)))?;
        let bounds: Vec<f64> = r.bounds.iter().flat_map(|b| b.iter().copied()).collect();
        out.push(list!(tile_index = r.tile_index, counts = r.counts, bounds = bounds).into());
    }
    Ok(List::from_values(out))
}

extendr_module! {
    mod ateratxrs;
    fn atera_tx_cells_rs;
}
