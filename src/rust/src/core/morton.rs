//! GiottoDisk's per-tile row order (`.morton_encode()` + `order()`).

/// `.morton_encode()`, with R's signed 32-bit wrap; `None` is NA (a
/// zero-width extent gives 0 / 0).
pub fn morton(x: f64, y: f64, e: [f64; 4]) -> Option<i32> {
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

/// Extent of the rows (xmin, xmax, ymin, ymax).
pub fn extent(x: &[f64], y: &[f64]) -> [f64; 4] {
    let fold = |v: &[f64], f: fn(f64, f64) -> f64, i: f64| v.iter().fold(i, |p, c| f(p, *c));
    [
        fold(x, f64::min, f64::INFINITY), fold(x, f64::max, f64::NEG_INFINITY),
        fold(y, f64::min, f64::INFINITY), fold(y, f64::max, f64::NEG_INFINITY),
    ]
}

/// Row order of a tile: Morton keys over the tile's own extent, stable, NA
/// last (R's `order()`). Returns the order and the extent.
pub fn order(x: &[f64], y: &[f64]) -> (Vec<u32>, [f64; 4]) {
    let e = extent(x, y);
    let keys: Vec<Option<i32>> = x.iter().zip(y).map(|(a, b)| morton(*a, *b, e)).collect();
    let mut ord: Vec<u32> = (0..x.len() as u32).collect();
    ord.sort_by_key(|&i| (keys[i as usize].is_none(), keys[i as usize].unwrap_or(0)));
    (ord, e)
}
