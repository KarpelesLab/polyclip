//! Measurements and point/region queries.

use crate::geom::Point;

/// Twice the signed area of a ring given as a vertex slice (exact, shoelace formula).
/// Positive for counter-clockwise rings.
pub fn ring_area2(pts: &[Point]) -> i128 {
    let n = pts.len();
    if n < 3 {
        return 0;
    }
    let o = pts[0];
    let mut s: i128 = 0;
    for i in 1..n - 1 {
        s += crate::predicates::orient(o, pts[i], pts[i + 1]);
    }
    s
}
