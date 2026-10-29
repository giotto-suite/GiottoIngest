use extendr_api::prelude::*;
use extendr_api::Result;

mod points;
mod polygons;
mod source;
mod tiles;

use crate::points::WriteOpts;
use crate::source::{Filter, Spec};

fn strings(s: &Strings) -> Vec<String> {
    s.iter().map(|x| x.to_string()).collect()
}

fn reals(r: &Robj, what: &str) -> Result<Vec<f64>> {
    r.as_real_slice().map(|v| v.to_vec())
        .ok_or_else(|| Error::Other(format!("`{what}` must be a double vector")))
}

fn spec(files: &Strings, x: &str, y: &str, flip_y: bool, f_cols: &Strings, f_ops: &Strings,
        f_vals: &Robj) -> Result<Spec> {
    let vals = reals(f_vals, "f_vals")?;
    let filters = f_cols.iter().zip(f_ops.iter()).zip(vals)
        .map(|((c, o), v)| Filter { col: c.to_string(), op: o.to_string(), value: v })
        .collect();
    Ok(Spec { files: strings(files), x: x.into(), y: y.into(), flip_y, filters })
}

fn bounds(b: &Robj) -> Result<Vec<[f64; 4]>> {
    let v = reals(b, "bounds")?;
    Ok(v.chunks(4).map(|c| [c[0], c[1], c[2], c[3]]).collect())
}

fn pool<T: Send>(threads: i32, f: impl FnOnce() -> std::result::Result<T, String> + Send) -> Result<T> {
    let p = rayon::ThreadPoolBuilder::new().num_threads(threads.max(1) as usize).build()
        .map_err(|e| Error::Other(e.to_string()))?;
    p.install(f).map_err(|e| Error::Other(format!("[giottoingest] {e}")))
}

/// Rows and extent of the filtered source.
/// @keywords internal
#[extendr]
fn gi_extent_rs(files: Strings, x: &str, y: &str, flip_y: bool, f_cols: Strings, f_ops: Strings,
                f_vals: Robj, threads: i32) -> Result<List> {
    let s = spec(&files, x, y, flip_y, &f_cols, &f_ops, &f_vals)?;
    let e = pool(threads, || points::extent(&s, threads as usize))?;
    Ok(list!(n = e.n, ext = e.e.to_vec()))
}

/// Rows per tile; `bounds` row-major xmin, xmax, ymin, ymax.
/// @keywords internal
#[extendr]
fn gi_count_tiles_rs(files: Strings, x: &str, y: &str, flip_y: bool, f_cols: Strings,
                     f_ops: Strings, f_vals: Robj, bounds_v: Robj, threads: i32) -> Result<Vec<f64>> {
    let s = spec(&files, x, y, flip_y, &f_cols, &f_ops, &f_vals)?;
    let b = bounds(&bounds_v)?;
    pool(threads, || points::count_tiles(&s, b, threads as usize))
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
    pool(threads, || points::write_tiles(&s, b, &o, threads as usize))
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

extendr_module! {
    mod giottoingest;
    fn gi_route_polygons_rs;
    fn gi_extent_rs;
    fn gi_count_tiles_rs;
    fn gi_write_point_tiles_rs;
}
