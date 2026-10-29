//! Tile membership exactly as GiottoDisk's contiguous `getTile()` on a
//! `freeTilePlan` filters: left and top edges inclusive; bottom inclusive
//! only for tiles on the set's lowest ymin, right only for tiles on its
//! highest xmax. The "set" is whatever bounds are passed, as in R, where
//! each `tileApply` call computes these limits over its own plan.

pub struct Tiles {
    b: Vec<[f64; 4]>, // xmin, xmax, ymin, ymax
    incl_bottom: Vec<bool>,
    incl_right: Vec<bool>,
    // grid index over the union of the tiles
    x0: f64,
    y0: f64,
    cw: f64,
    ch: f64,
    nx: usize,
    ny: usize,
    cells: Vec<Vec<u32>>,
}

impl Tiles {
    /// `x_f32` / `y_f32`: the coordinate column is float32. Arrow then casts
    /// the bound to float32 before comparing, so the bounds are rounded the
    /// same way here; the coordinates themselves are exact f32 values.
    pub fn new(b: Vec<[f64; 4]>, x_f32: bool, y_f32: bool) -> Tiles {
        // which edges are inclusive: decided in R on the double bounds
        let gx = b.iter().map(|t| t[1]).fold(f64::NEG_INFINITY, f64::max);
        let gy = b.iter().map(|t| t[2]).fold(f64::INFINITY, f64::min);
        let incl_bottom = b.iter().map(|t| t[2] == gy).collect();
        let incl_right = b.iter().map(|t| t[1] == gx).collect();
        // the comparisons themselves: bounds as the column type sees them
        let r = |v: f64, on: bool| if on { v as f32 as f64 } else { v };
        let b: Vec<[f64; 4]> = b.iter()
            .map(|t| [r(t[0], x_f32), r(t[1], x_f32), r(t[2], y_f32), r(t[3], y_f32)])
            .collect();
        let gx = b.iter().map(|t| t[1]).fold(f64::NEG_INFINITY, f64::max);
        let gy = b.iter().map(|t| t[2]).fold(f64::INFINITY, f64::min);
        let x0 = b.iter().map(|t| t[0]).fold(f64::INFINITY, f64::min);
        let y1 = b.iter().map(|t| t[3]).fold(f64::NEG_INFINITY, f64::max);
        let n = ((b.len() as f64).sqrt().ceil() as usize * 2).clamp(1, 1024);
        let (nx, ny) = (n, n);
        let cw = ((gx - x0) / nx as f64).max(f64::MIN_POSITIVE);
        let ch = ((y1 - gy) / ny as f64).max(f64::MIN_POSITIVE);
        let mut cells = vec![Vec::new(); nx * ny];
        let cell = |v: f64, o: f64, w: f64, n: usize| -> usize {
            let k = ((v - o) / w).floor();
            if k < 0.0 { 0 } else { (k as usize).min(n - 1) }
        };
        for (k, t) in b.iter().enumerate() {
            let (i0, i1) = (cell(t[0], x0, cw, nx), cell(t[1], x0, cw, nx));
            let (j0, j1) = (cell(t[2], gy, ch, ny), cell(t[3], gy, ch, ny));
            for j in j0..=j1 {
                for i in i0..=i1 {
                    cells[j * nx + i].push(k as u32);
                }
            }
        }
        Tiles { b, incl_bottom, incl_right, x0, y0: gy, cw, ch, nx, ny, cells }
    }

    pub fn len(&self) -> usize { self.b.len() }

    fn contains(&self, k: usize, x: f64, y: f64) -> bool {
        let t = &self.b[k];
        let bottom = if self.incl_bottom[k] { y >= t[2] } else { y > t[2] };
        let right = if self.incl_right[k] { x <= t[1] } else { x < t[1] };
        bottom && x >= t[0] && y <= t[3] && right
    }

    /// Every tile containing (x, y), in tile order. R filters each tile
    /// independently, so a point that met two tiles' filters would be in both.
    pub fn matches(&self, x: f64, y: f64, out: &mut Vec<u32>) {
        out.clear();
        if !(x.is_finite() && y.is_finite()) {
            return;
        }
        let ci = ((x - self.x0) / self.cw).floor();
        let cj = ((y - self.y0) / self.ch).floor();
        // a point on the far edge falls one cell out; clamp, then the exact
        // test decides
        let i = if ci < 0.0 { 0 } else { (ci as usize).min(self.nx - 1) };
        let j = if cj < 0.0 { 0 } else { (cj as usize).min(self.ny - 1) };
        for &k in &self.cells[j * self.nx + i] {
            if self.contains(k as usize, x, y) {
                out.push(k);
            }
        }
        out.sort_unstable();
        out.dedup();
    }
}
