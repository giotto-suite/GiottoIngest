use extendr_api::prelude::*;
use extendr_api::Result;

// The system allocator on macOS turns the writers' large reallocations into
// kernel `vm_copy` calls that contend across threads (a 2 s Xenium write took
// 3-360 s); mimalloc keeps them in user space.
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

mod core;
mod input;
mod points;
mod polygons;
mod tiles;

use crate::input::parquet::{Filter, Spec};
use crate::input::zarr::atera::{run_cell, CellPlan, CellResult, Lookups, Opts};
use crate::points::WriteOpts;

fn strings(s: &Strings) -> Vec<String> {
    s.iter().map(|x| x.to_string()).collect()
}

fn reals(r: &Robj, what: &str) -> Result<Vec<f64>> {
    r.as_real_slice().map(|v| v.to_vec())
        .ok_or_else(|| Error::Other(format!("`{what}` must be a double vector")))
}

/// `f_vals` is a list with one double vector per filter (several for `in`).
fn spec(files: &Strings, x: &str, y: &str, flip_y: bool, f_cols: &Strings, f_ops: &Strings,
        f_vals: &Robj) -> Result<Spec> {
    let vals = f_vals.as_list().ok_or_else(|| Error::Other("`f_vals` must be a list".into()))?;
    let mut filters = Vec::with_capacity(f_cols.len());
    for (k, (c, o)) in f_cols.iter().zip(f_ops.iter()).enumerate() {
        let v = reals(&vals.elt(k)?, "filter value")?;
        filters.push(Filter { col: c.to_string(), op: o.to_string(), values: v });
    }
    Ok(Spec { files: strings(files), x: x.into(), y: y.into(), flip_y, filters })
}

/// Threads to use: at least one (the R side rejects NA).
fn nthreads(threads: i32) -> usize {
    threads.max(1) as usize
}

fn bounds(b: &Robj) -> Result<Vec<[f64; 4]>> {
    let v = reals(b, "bounds")?;
    Ok(v.chunks(4).map(|c| [c[0], c[1], c[2], c[3]]).collect())
}

/// Runs `f` on a pool of `threads`, interruptible from R.
fn pool<T: Send>(threads: i32, f: impl FnOnce() -> std::result::Result<T, String> + Send) -> Result<T> {
    let p = rayon::ThreadPoolBuilder::new().num_threads(nthreads(threads)).build()
        .map_err(|e| Error::Other(e.to_string()))?;
    crate::core::cancel::run(|| p.install(f)).map_err(|e| Error::Other(format!("[GiottoIngest] {e}")))
}

/// Rows and extent of the filtered source.
/// @keywords internal
#[extendr]
fn gi_extent_rs(files: Strings, x: &str, y: &str, flip_y: bool, f_cols: Strings, f_ops: Strings,
                f_vals: Robj, threads: i32) -> Result<List> {
    let s = spec(&files, x, y, flip_y, &f_cols, &f_ops, &f_vals)?;
    let e = pool(threads, || points::extent(&s, nthreads(threads)))?;
    Ok(list!(n = e.n, ext = e.e.to_vec()))
}

/// Rows per tile; `bounds` row-major xmin, xmax, ymin, ymax.
/// @keywords internal
#[extendr]
fn gi_count_tiles_rs(files: Strings, x: &str, y: &str, flip_y: bool, f_cols: Strings,
                     f_ops: Strings, f_vals: Robj, bounds_v: Robj, threads: i32) -> Result<Vec<f64>> {
    let s = spec(&files, x, y, flip_y, &f_cols, &f_ops, &f_vals)?;
    let b = bounds(&bounds_v)?;
    pool(threads, || points::count_tiles(&s, b, nthreads(threads)))
}

/// Routes rows to tiles and writes one parquet file per tile.
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn gi_write_point_tiles_rs(files: Strings, x: &str, y: &str, flip_y: bool, f_cols: Strings,
                           f_ops: Strings, f_vals: Robj, bounds_v: Robj, tile_dirs: Strings,
                           id_col: &str, to_string: Strings, drop_cols: Strings, root: &str,
                           geo_json: &str, spill_dir: &str, spill_rows: f64,
                           geom: bool, threads: i32) -> Result<Vec<f64>> {
    let s = spec(&files, x, y, flip_y, &f_cols, &f_ops, &f_vals)?;
    let b = bounds(&bounds_v)?;
    let o = WriteOpts {
        id_col: id_col.into(), to_string: strings(&to_string), drop: strings(&drop_cols),
        root: root.into(), tile_dirs: strings(&tile_dirs), geo_json: geo_json.into(),
        spill_dir: spill_dir.into(), spill_rows: spill_rows as usize, geom,
    };
    pool(threads, || points::write_tiles(&s, b, &o, nthreads(threads)))
}

/// What this build can do, for GiottoDisk to gate on (see `has_capability()`).
/// Names are only ever added; a build without one falls back to R.
/// @keywords internal
#[extendr]
fn capabilities_rs() -> Vec<String> {
    ["points_parquet", "points_atera_zarr"].iter().map(|s| s.to_string()).collect()
}

/// GeoParquet metadata of a point store with extent `ext` (xmin, xmax, ymin,
/// ymax), as GiottoDisk's `.geoparquet_metadata()` writes it.
/// @keywords internal
#[extendr]
fn gi_geo_json_rs(ext: &[f64]) -> Result<String> {
    if ext.len() != 4 {
        return Err(Error::Other("`ext` must be xmin, xmax, ymin, ymax".into()));
    }
    Ok(crate::core::geo::points_geo([ext[0], ext[1], ext[2], ext[3]]))
}

/// Assigns polygons to tiles by envelope centroid and writes each tile's
/// vertex rows, in row_index order, to `out_dir/<tile_files>`.
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn gi_route_polygons_rs(files: Strings, x: &str, y: &str, flip_y: bool, f_cols: Strings,
                        f_ops: Strings, f_vals: Robj, bounds_v: Robj, group: &str,
                        out_dir: &str, tile_files: Strings, threads: i32) -> Result<List> {
    let s = spec(&files, x, y, flip_y, &f_cols, &f_ops, &f_vals)?;
    let b = bounds(&bounds_v)?;
    let tf = strings(&tile_files);
    let r = pool(threads, || polygons::route_polygons(&s, b, group, out_dir, &tf))?;
    Ok(list!(counts = r.counts, groups = r.groups, unassigned = r.unassigned))
}

/// Stream Atera grid cells into tile partitions on `threads` threads.
///
/// Internal; called by [gi_atera_tx_cells()], which validates the arguments.
/// Each thread takes the next grid cell in the given order, so pass them
/// largest first.
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
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
        gene_names: strings(&gene_names),
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
    crate::core::cancel::run(|| {
        std::thread::scope(|sc| {
            for _ in 0..nthreads(threads) {
                sc.spawn(|| loop {
                    let k = next.fetch_add(1, Ordering::SeqCst);
                    if k >= plans.len() || crate::core::cancel::check().is_err() {
                        break;
                    }
                    let r = run_cell(zarr_path, &plans[k], &lk, &opts);
                    results.lock().unwrap()[k] = Some(r);
                });
            }
        });
        Ok(())
    }).map_err(|e| Error::Other(format!("[GiottoIngest] {e}")))?;

    let mut out: Vec<Robj> = Vec::with_capacity(plans.len());
    for (k, r) in results.into_inner().unwrap().into_iter().enumerate() {
        let r = r.ok_or_else(|| Error::Other(format!("cell {} not run", plans[k].tile)))?
            .map_err(|e| Error::Other(format!("[GiottoIngest] Atera cell {}: {e}", plans[k].tile)))?;
        let bounds: Vec<f64> = r.bounds.iter().flat_map(|b| b.iter().copied()).collect();
        out.push(list!(tile_index = r.tile_index, counts = r.counts, bounds = bounds).into());
    }
    Ok(List::from_values(out))
}

extendr_module! {
    mod GiottoIngest;
    fn atera_tx_cells_rs;
    fn capabilities_rs;
    fn gi_geo_json_rs;
    fn gi_route_polygons_rs;
    fn gi_extent_rs;
    fn gi_count_tiles_rs;
    fn gi_write_point_tiles_rs;
}
