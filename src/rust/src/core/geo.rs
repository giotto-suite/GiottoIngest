//! GeoParquet `geo` metadata exactly as GiottoDisk's
//! `.geoparquet_metadata(geom_col = "geom", geomtype = "points", crs = "")`
//! writes it: `jsonlite::toJSON(auto_unbox = TRUE)` with its default
//! `digits = 4` number formatting.

/// A number as `jsonlite::toJSON(x, digits = 4)` writes it.
pub fn json_num(x: f64) -> String {
    if x.is_nan() {
        return "\"NA\"".into(); // never in a bbox; jsonlite's string form
    }
    let neg = x.is_sign_negative();
    let a = x.abs();
    let body = if a == 0.0 {
        "0".to_string()
    } else if a <= 1e-5 || a > 1e16 {
        sci(a)
    } else {
        fixed4(a)
    };
    if neg { format!("-{body}") } else { body }
}

/// C's `%.4e` with trailing mantissa zeros dropped: `1e-05`, `6.6667e-06`.
fn sci(a: f64) -> String {
    let s = format!("{a:.4e}");
    let (m, e) = s.split_once('e').unwrap();
    let m = if m.contains('.') { m.trim_end_matches('0').trim_end_matches('.') } else { m };
    let e: i32 = e.parse().unwrap();
    format!("{m}e{}{:02}", if e < 0 { '-' } else { '+' }, e.abs())
}

/// Four decimals, round half to even on the scaled fraction, zeros trimmed.
fn fixed4(a: f64) -> String {
    let mut whole = a.trunc();
    let frac = (a - whole) * 10_000.0;
    let mut f = frac.trunc() as u64;
    let diff = frac - f as f64;
    if diff > 0.5 || (diff == 0.5 && f & 1 == 1) {
        f += 1;
    }
    if f >= 10_000 {
        f = 0;
        whole += 1.0;
    }
    let w = format!("{whole:.0}");
    if f == 0 {
        w
    } else {
        let d = format!("{f:04}");
        format!("{w}.{}", d.trim_end_matches('0'))
    }
}

/// The point store's `geo` metadata for an extent (xmin, xmax, ymin, ymax).
pub fn points_geo(e: [f64; 4]) -> String {
    format!(
        "{{\"version\":\"1.0.0\",\"primary_column\":\"geom\",\"columns\":{{\"geom\":{{\"encoding\":\"WKB\",\"geometry_types\":[\"Point\"],\"crs\":null,\"edges\":\"planar\",\"bbox\":[{},{},{},{}]}}}}}}",
        json_num(e[0]), json_num(e[2]), json_num(e[1]), json_num(e[3])
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    // values and jsonlite::toJSON(x, digits = 4) output, probed in R
    #[test]
    fn matches_jsonlite() {
        let cases: &[(f64, &str)] = &[
            (0.0, "0"), (-0.0, "-0"), (1.0, "1"), (0.1, "0.1"), (0.00015, "0.0001"),
            (0.0001, "0.0001"), (9.999e-05, "0.0001"), (5e-05, "0"), (1.234567e-05, "0"),
            (1.5e-05, "0"), (1e-05, "1e-05"), (9.99999e-06, "1e-05"), (1e-15, "1e-15"),
            (123456789.123456, "123456789.1235"), (1e15, "1000000000000000"),
            (1.5e15, "1500000000000000"), (2147483648.0, "2147483648"),
            (2147483647.5, "2147483647.5"), (1e21, "1e+21"), (-2.5e-06, "-2.5e-06"),
            (0.99999, "1"), (0.000123456, "0.0001"), (12.3, "12.3"),
            (5e-324, "4.9407e-324"), (1.0 / 3.0, "0.3333"), (2.0 / 3.0 * 1e-5, "6.6667e-06"),
            (1.00000001e-5, "0"), (1.1e-5, "0"), (0.00025, "0.0002"), (0.00035, "0.0004"),
            (0.00045, "0.0004"), (1.00005, "1.0001"), (2.00015, "2.0002"), (0.5, "0.5"),
            (1e16, "10000000000000000"), (1e17, "1e+17"), (9e20, "9e+20"), (1e20, "1e+20"),
            (123456789012.12345, "123456789012.1234"), (4294967296.5, "4294967296.5"),
            (4294967295.99995, "4294967295.9999"), (1e9 + 0.00005, "1000000000"),
            (-0.00025, "-0.0002"), (-1.00005, "-1.0001"), (1e-300, "1e-300"),
            (1.23456789, "1.2346"), (-3.00001, "-3"), (4e6, "4000000"),
            (40598.387123, "40598.3871"), (-22706.3312, "-22706.3312"), (2063.89, "2063.89"),
        ];
        for (x, want) in cases {
            assert_eq!(json_num(*x), *want, "x = {x:e}");
        }
    }

    #[test]
    fn geo_string() {
        assert_eq!(
            points_geo([1.23456789, 2.0, -3.00001, 4e6]),
            "{\"version\":\"1.0.0\",\"primary_column\":\"geom\",\"columns\":{\"geom\":{\"encoding\":\"WKB\",\"geometry_types\":[\"Point\"],\"crs\":null,\"edges\":\"planar\",\"bbox\":[1.2346,-3,2,4000000]}}}"
        );
    }
}
