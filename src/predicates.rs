//! Exact geometric predicates on integer points.
//!
//! All functions assume coordinates within `±MAX_COORD` (2^40); every intermediate value
//! then fits in `i128`.

use crate::geom::{Point, Rect};
use core::cmp::Ordering;

/// Cross product of `b - a` and `c - a` (twice the signed area of triangle `abc`).
/// Positive when `c` is to the left of the directed line `a -> b`.
#[inline]
pub fn orient(a: Point, b: Point, c: Point) -> i128 {
    let abx = (b.x - a.x) as i128;
    let aby = (b.y - a.y) as i128;
    let acx = (c.x - a.x) as i128;
    let acy = (c.y - a.y) as i128;
    abx * acy - aby * acx
}

/// Cross product of two vectors given as points.
#[inline]
pub fn cross(u: Point, v: Point) -> i128 {
    u.x as i128 * v.y as i128 - u.y as i128 * v.x as i128
}

/// Dot product of two vectors given as points.
#[inline]
pub fn dot(u: Point, v: Point) -> i128 {
    u.x as i128 * v.x as i128 + u.y as i128 * v.y as i128
}

/// `b - a` as a vector.
#[inline]
pub fn sub(b: Point, a: Point) -> Point {
    Point::new(b.x - a.x, b.y - a.y)
}

/// Squared length of `b - a`.
#[inline]
pub fn dist2(a: Point, b: Point) -> i128 {
    let dx = (b.x - a.x) as i128;
    let dy = (b.y - a.y) as i128;
    dx * dx + dy * dy
}

/// `true` when `p` lies on the closed segment `a-b`.
#[inline]
pub fn on_segment(a: Point, b: Point, p: Point) -> bool {
    orient(a, b, p) == 0
        && p.x >= a.x.min(b.x)
        && p.x <= a.x.max(b.x)
        && p.y >= a.y.min(b.y)
        && p.y <= a.y.max(b.y)
}

/// `true` when `p` lies on the segment `a-b` but is neither endpoint.
#[inline]
pub fn in_segment_interior(a: Point, b: Point, p: Point) -> bool {
    p != a && p != b && on_segment(a, b, p)
}

/// `true` when the closed segments `a-b` and `c-d` share at least one point.
pub fn segments_intersect(a: Point, b: Point, c: Point, d: Point) -> bool {
    if a.x.max(b.x) < c.x.min(d.x)
        || c.x.max(d.x) < a.x.min(b.x)
        || a.y.max(b.y) < c.y.min(d.y)
        || c.y.max(d.y) < a.y.min(b.y)
    {
        return false;
    }
    let o1 = orient(a, b, c).signum();
    let o2 = orient(a, b, d).signum();
    let o3 = orient(c, d, a).signum();
    let o4 = orient(c, d, b).signum();
    if o1 * o2 < 0 && o3 * o4 < 0 {
        return true;
    }
    (o1 == 0 && on_segment(a, b, c))
        || (o2 == 0 && on_segment(a, b, d))
        || (o3 == 0 && on_segment(c, d, a))
        || (o4 == 0 && on_segment(c, d, b))
}

/// `true` when the segments cross at a single point interior to both (a transversal
/// crossing, no endpoint involved).
#[inline]
pub fn segments_cross_properly(a: Point, b: Point, c: Point, d: Point) -> bool {
    let o1 = orient(a, b, c).signum();
    let o2 = orient(a, b, d).signum();
    if o1 * o2 >= 0 {
        return false;
    }
    let o3 = orient(c, d, a).signum();
    let o4 = orient(c, d, b).signum();
    o3 * o4 < 0
}

/// Floor division for `i128` with positive divisor.
#[inline]
pub fn floor_div(a: i128, b: i128) -> i128 {
    debug_assert!(b > 0);
    let q = a / b;
    if (a % b) < 0 { q - 1 } else { q }
}

/// For two properly crossing segments, the crossing point rounded to the nearest integer
/// point (ties rounded towards +infinity, i.e. `floor(v + 1/2)`), which is the center of the
/// half-open unit pixel containing the exact crossing.
pub fn rounded_crossing(a: Point, b: Point, c: Point, d: Point) -> Point {
    // a + (b - a) * t,  t = o3 / (o3 - o4) where o3 = orient(c,d,a), o4 = orient(c,d,b).
    let o3 = orient(c, d, a);
    let o4 = orient(c, d, b);
    let mut num = o3;
    let mut den = o3 - o4;
    if den < 0 {
        num = -num;
        den = -den;
    }
    let dx = (b.x - a.x) as i128;
    let dy = (b.y - a.y) as i128;
    let rx = floor_div(2 * num * dx + den, 2 * den);
    let ry = floor_div(2 * num * dy + den, 2 * den);
    Point::new(a.x + rx as i64, a.y + ry as i64)
}

/// The exact crossing point of two non-parallel lines as `f64` (for diagnostics only).
pub fn crossing_f64(a: Point, b: Point, c: Point, d: Point) -> (f64, f64) {
    let o3 = orient(c, d, a);
    let o4 = orient(c, d, b);
    let den = o3 - o4;
    if den == 0 {
        return (a.x as f64, a.y as f64);
    }
    let t = o3 as f64 / den as f64;
    (
        a.x as f64 + (b.x - a.x) as f64 * t,
        a.y as f64 + (b.y - a.y) as f64 * t,
    )
}

/// A non-negative rational number `n / d` with `d > 0`, compared exactly.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Frac {
    pub n: i128,
    pub d: i128,
}

impl Frac {
    #[inline]
    pub fn cmp(self, o: Frac) -> Ordering {
        (self.n * o.d).cmp(&(o.n * self.d))
    }
}

/// One side of a parametric interval: value and whether the bound is excluded.
#[derive(Clone, Copy, Debug)]
struct Bound {
    v: Frac,
    open: bool,
}

/// Clips the parametric interval of segment `p + t * dv` (in `[lo, hi]`) against
/// `min <= coord < max` (or `<= max` when `max_closed`). Returns `false` when empty.
#[inline]
fn clip_axis(
    p: i128,
    dv: i128,
    min: i128,
    max: i128,
    max_closed: bool,
    lo: &mut Bound,
    hi: &mut Bound,
) -> bool {
    if dv == 0 {
        return p >= min && (if max_closed { p <= max } else { p < max });
    }
    // Constraint min <= p + t dv : closed. Constraint p + t dv < max (or <=).
    let (lower, upper) = if dv > 0 {
        (
            Bound {
                v: Frac { n: min - p, d: dv },
                open: false,
            },
            Bound {
                v: Frac { n: max - p, d: dv },
                open: !max_closed,
            },
        )
    } else {
        // t <= (p - min) / -dv (closed);  t > (p - max) / -dv (open unless closed).
        (
            Bound {
                v: Frac { n: p - max, d: -dv },
                open: !max_closed,
            },
            Bound {
                v: Frac { n: p - min, d: -dv },
                open: false,
            },
        )
    };
    match lower.v.cmp(lo.v) {
        Ordering::Greater => *lo = lower,
        Ordering::Equal => lo.open |= lower.open,
        Ordering::Less => {}
    }
    match upper.v.cmp(hi.v) {
        Ordering::Less => *hi = upper,
        Ordering::Equal => hi.open |= upper.open,
        Ordering::Greater => {}
    }
    match lo.v.cmp(hi.v) {
        Ordering::Less => true,
        Ordering::Equal => !lo.open && !hi.open,
        Ordering::Greater => false,
    }
}

/// If segment `a-b` meets the half-open pixel `[c.x - 1/2, c.x + 1/2) x [c.y - 1/2, c.y + 1/2)`,
/// returns the parameter at which it enters the pixel (a fraction of the segment, as
/// `(num, den)` with `den > 0`) and whether that bound is open (the segment touches the
/// entry value only from above).
pub(crate) fn segment_pixel_entry(a: Point, b: Point, c: Point) -> Option<(Frac, bool)> {
    // Work in doubled coordinates so that pixel bounds are integers.
    let px = 2 * a.x as i128;
    let py = 2 * a.y as i128;
    let dx = 2 * (b.x - a.x) as i128;
    let dy = 2 * (b.y - a.y) as i128;
    let mut lo = Bound {
        v: Frac { n: 0, d: 1 },
        open: false,
    };
    let mut hi = Bound {
        v: Frac { n: 1, d: 1 },
        open: false,
    };
    let cx = 2 * c.x as i128;
    let cy = 2 * c.y as i128;
    if !clip_axis(px, dx, cx - 1, cx + 1, false, &mut lo, &mut hi) {
        return None;
    }
    if !clip_axis(py, dy, cy - 1, cy + 1, false, &mut lo, &mut hi) {
        return None;
    }
    Some((lo.v, lo.open))
}

/// `true` when segment `a-b` meets the closed rectangle `r`.
pub fn segment_meets_rect(a: Point, b: Point, r: &Rect) -> bool {
    if a.x.max(b.x) < r.min.x
        || a.x.min(b.x) > r.max.x
        || a.y.max(b.y) < r.min.y
        || a.y.min(b.y) > r.max.y
    {
        return false;
    }
    if r.contains_point(a) || r.contains_point(b) {
        return true;
    }
    // Bounding boxes overlap: the segment meets the rectangle unless all four corners are
    // strictly on the same side of its supporting line.
    let corners = [
        r.min,
        Point::new(r.max.x, r.min.y),
        r.max,
        Point::new(r.min.x, r.max.y),
    ];
    let mut pos = false;
    let mut neg = false;
    for c in corners {
        let o = orient(a, b, c);
        if o == 0 {
            return true;
        }
        if o > 0 {
            pos = true;
        } else {
            neg = true;
        }
    }
    pos && neg
}

/// Angular comparison of two direction vectors that both point into the half-plane
/// `x > 0 || (x == 0 && y > 0)`: orders them bottom to top (counter-clockwise).
#[inline]
pub(crate) fn cmp_dir_halfplane(u: Point, v: Point) -> Ordering {
    0.cmp(&cross(u, v))
}

/// Full-circle angle comparison of non-zero direction vectors, counter-clockwise starting
/// from the positive x axis (angle 0 inclusive).
pub(crate) fn cmp_angle(u: Point, v: Point) -> Ordering {
    #[inline]
    fn half(p: Point) -> u8 {
        // 0 for angles in [0, pi), 1 for [pi, 2pi).
        if p.y > 0 || (p.y == 0 && p.x > 0) {
            0
        } else {
            1
        }
    }
    half(u).cmp(&half(v)).then_with(|| 0.cmp(&cross(u, v)))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(x: i64, y: i64) -> Point {
        Point::new(x, y)
    }

    #[test]
    fn orientation() {
        assert!(orient(p(0, 0), p(10, 0), p(5, 1)) > 0);
        assert!(orient(p(0, 0), p(10, 0), p(5, -1)) < 0);
        assert_eq!(orient(p(0, 0), p(10, 0), p(20, 0)), 0);
    }

    #[test]
    fn crossing_rounding() {
        // Lines y = x and y = -x + 3 cross at (1.5, 1.5) -> rounds to (2, 2).
        assert_eq!(
            rounded_crossing(p(0, 0), p(4, 4), p(0, 3), p(3, 0)),
            p(2, 2)
        );
        // y = x and y = -x + 1: (0.5,0.5) -> (1,1).
        assert_eq!(
            rounded_crossing(p(-2, -2), p(2, 2), p(-1, 2), p(2, -1)),
            p(1, 1)
        );
        // Exact at extreme range.
        let m = crate::MAX_COORD;
        let c = rounded_crossing(p(-m, -m), p(m, m), p(-m, m), p(m, -m));
        assert_eq!(c, p(0, 0));
    }

    #[test]
    fn pixel_entry() {
        // Horizontal segment along y = 0 through pixels (0,0) .. (4,0).
        assert!(segment_pixel_entry(p(0, 0), p(4, 0), p(2, 0)).is_some());
        assert!(segment_pixel_entry(p(0, 0), p(4, 0), p(2, 1)).is_none());
        // Segment along y = x + 1/2 would graze; use doubled-safe example: (0,0)-(2,1) passes
        // through (1, 0.5) which belongs to pixel (1,1) (half-open, upper side closed below).
        assert!(segment_pixel_entry(p(0, 0), p(2, 1), p(1, 1)).is_some());
        assert!(segment_pixel_entry(p(0, 0), p(2, 1), p(1, 0)).is_some());
        // Segment (0,0)-(1,1) touches the corner (0.5,0.5) which belongs to pixel (1,1) only.
        assert!(segment_pixel_entry(p(0, 0), p(1, 1), p(1, 0)).is_none());
        assert!(segment_pixel_entry(p(0, 0), p(1, 1), p(0, 1)).is_none());
    }

    #[test]
    fn rect_meet() {
        let r = Rect::new(p(0, 0), p(10, 10));
        assert!(segment_meets_rect(p(-5, 5), p(15, 5), &r));
        assert!(segment_meets_rect(p(-5, -5), p(0, 0), &r));
        assert!(!segment_meets_rect(p(-5, 6), p(6, 17), &r));
        assert!(segment_meets_rect(p(-5, 5), p(5, 15), &r));
    }

    #[test]
    fn angles() {
        let dirs = [
            p(1, 0),
            p(1, 1),
            p(0, 1),
            p(-1, 1),
            p(-1, 0),
            p(-1, -1),
            p(0, -1),
            p(1, -1),
        ];
        for i in 0..dirs.len() {
            for j in 0..dirs.len() {
                assert_eq!(cmp_angle(dirs[i], dirs[j]), i.cmp(&j), "{i} {j}");
            }
        }
    }
}
