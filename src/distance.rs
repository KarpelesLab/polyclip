//! Exact minimum distance between geometries.
//!
//! Distances are between point sets: areal geometries are closed regions, so a point inside
//! a polygon is at distance zero from it, and two overlapping polygons are at distance zero.
//! Squared distances are exact rationals ([`SqDist`]); comparisons against a threshold are
//! exact as well.

use crate::geom::{Point, PointF, Rect};
use crate::predicates::{cross, crossing_f64, dist2, dot, on_segment, segments_intersect, sub};
use crate::query::{Geometry, any_pair, collect_segments};
use crate::wide::{U384, cmp_products};
use core::cmp::Ordering;

/// An exact squared distance: the rational `num / den`.
#[derive(Clone, Copy, Debug)]
pub struct SqDist {
    num: U384,
    den: u128,
}

impl SqDist {
    /// Zero.
    pub const ZERO: SqDist = SqDist {
        num: U384([0; 6]),
        den: 1,
    };

    fn int(v: u128) -> SqDist {
        SqDist {
            num: U384::from_u128(v),
            den: 1,
        }
    }

    /// `c^2 / den`.
    fn ratio_sq(c: u128, den: u128) -> SqDist {
        SqDist {
            num: U384::mul_u128(c, c),
            den,
        }
    }

    /// `true` when the distance is zero.
    pub fn is_zero(&self) -> bool {
        self.num == U384::default()
    }

    /// The squared distance as `f64` (correctly ordered up to rounding).
    pub fn to_f64(&self) -> f64 {
        self.num.to_f64() / self.den as f64
    }

    /// The distance as `f64`.
    pub fn distance_f64(&self) -> f64 {
        libm::sqrt(self.to_f64())
    }

    /// Exact comparison of this squared distance with `d * d`.
    pub fn cmp_dist(&self, d: u64) -> Ordering {
        let d2 = d as u128 * d as u128;
        // num / den vs d2  <=>  num vs d2 * den
        self.num.cmp(&U384::mul_u128(d2, self.den))
    }
}

impl PartialEq for SqDist {
    fn eq(&self, o: &Self) -> bool {
        self.cmp(o) == Ordering::Equal
    }
}

impl Eq for SqDist {}

impl PartialOrd for SqDist {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}

impl Ord for SqDist {
    fn cmp(&self, o: &Self) -> Ordering {
        if self.den == o.den {
            return self.num.cmp(&o.num);
        }
        self.num.mul_small(o.den).cmp(&o.num.mul_small(self.den))
    }
}

/// Result of [`distance`]: the exact squared distance and a pair of closest points (one on
/// each geometry, as `f64` since they need not be integer points).
#[derive(Clone, Copy, Debug)]
pub struct Closest {
    /// Exact squared distance.
    pub sq: SqDist,
    /// Closest point on the first geometry.
    pub a: PointF,
    /// Closest point on the second geometry.
    pub b: PointF,
}

/// Exact squared distance from `p` to segment `a-b`, and the closest point of the segment.
fn point_segment(p: Point, a: Point, b: Point) -> (SqDist, PointF) {
    if a == b {
        return (SqDist::int(dist2(p, a) as u128), a.into());
    }
    let d = sub(b, a);
    let t = dot(sub(p, a), d);
    if t <= 0 {
        return (SqDist::int(dist2(p, a) as u128), a.into());
    }
    let len2 = dot(d, d);
    if t >= len2 {
        return (SqDist::int(dist2(p, b) as u128), b.into());
    }
    let c = cross(d, sub(p, a)).unsigned_abs();
    let f = t as f64 / len2 as f64;
    let q = PointF::new(a.x as f64 + d.x as f64 * f, a.y as f64 + d.y as f64 * f);
    (SqDist::ratio_sq(c, len2 as u128), q)
}

/// `true` when the distance from `p` to segment `a-b` is below `d` (`d2 = d * d`).
#[inline]
fn point_segment_lt(p: Point, a: Point, b: Point, d2: u128) -> bool {
    if a == b {
        return (dist2(p, a) as u128) < d2;
    }
    let d = sub(b, a);
    let t = dot(sub(p, a), d);
    if t <= 0 {
        return (dist2(p, a) as u128) < d2;
    }
    let len2 = dot(d, d);
    if t >= len2 {
        return (dist2(p, b) as u128) < d2;
    }
    let c = cross(d, sub(p, a)).unsigned_abs();
    cmp_products(c, c, d2, len2 as u128) == Ordering::Less
}

/// A point common to two intersecting segments (for reporting).
fn common_point(a: Point, b: Point, c: Point, d: Point) -> PointF {
    for (p, (s, t)) in [(c, (a, b)), (d, (a, b)), (a, (c, d)), (b, (c, d))] {
        if on_segment(s, t, p) {
            return p.into();
        }
    }
    let (x, y) = crossing_f64(a, b, c, d);
    PointF::new(x, y)
}

/// Exact distance between two segments with closest points.
fn segment_segment(a: Point, b: Point, c: Point, d: Point) -> (SqDist, PointF, PointF) {
    if segments_intersect(a, b, c, d) {
        let p = common_point(a, b, c, d);
        return (SqDist::ZERO, p, p);
    }
    let mut best: Option<(SqDist, PointF, PointF)> = None;
    let mut consider = |s: SqDist, x: PointF, y: PointF| {
        if best.as_ref().is_none_or(|b| s < b.0) {
            best = Some((s, x, y));
        }
    };
    let (s, q) = point_segment(a, c, d);
    consider(s, a.into(), q);
    let (s, q) = point_segment(b, c, d);
    consider(s, b.into(), q);
    let (s, q) = point_segment(c, a, b);
    consider(s, q, c.into());
    let (s, q) = point_segment(d, a, b);
    consider(s, q, d.into());
    best.unwrap()
}

/// `true` when two segments are closer than `d` (`d2 = d * d`).
#[inline]
fn segment_segment_lt(a: Point, b: Point, c: Point, d: Point, d2: u128) -> bool {
    segments_intersect(a, b, c, d)
        || point_segment_lt(a, c, d, d2)
        || point_segment_lt(b, c, d, d2)
        || point_segment_lt(c, a, b, d2)
        || point_segment_lt(d, a, b, d2)
}

/// Squared gap between two rectangles (0 when they overlap).
fn rect_gap2(a: &Rect, b: &Rect) -> u128 {
    let dx = (b.min.x - a.max.x).max(a.min.x - b.max.x).max(0) as u128;
    let dy = (b.min.y - a.max.y).max(a.min.y - b.max.y).max(0) as u128;
    dx * dx + dy * dy
}

fn containment<A: Geometry + ?Sized, B: Geometry + ?Sized>(outer: &A, inner: &B) -> Option<Point> {
    crate::query::any_component_inside(outer, inner)
}

/// Exact minimum distance between two geometries, with a pair of closest points. `None`
/// when either geometry is empty or has coordinates outside `±`[`MAX_COORD`](crate::MAX_COORD).
///
/// ```
/// use polyclip::{distance, Point, Ring};
/// let a = Ring::from([(0, 0), (10, 0), (10, 10), (0, 10)]);
/// let b = Ring::from([(13, 14), (20, 14), (20, 20)]);
/// let c = distance(&a, &b).unwrap();
/// assert_eq!(c.sq.cmp_dist(5), core::cmp::Ordering::Equal); // (3, 4) apart
/// assert_eq!(distance(&a, &Point::new(5, 5)).unwrap().sq.is_zero(), true);
/// ```
pub fn distance<A: Geometry + ?Sized, B: Geometry + ?Sized>(a: &A, b: &B) -> Option<Closest> {
    let ba = a.bbox()?;
    let bb = b.bbox()?;
    if ![ba.min, ba.max, bb.min, bb.max]
        .iter()
        .all(|p| p.in_range())
    {
        return None;
    }
    let everything = Rect {
        min: Point::new(i64::MIN, i64::MIN),
        max: Point::new(i64::MAX, i64::MAX),
    };
    let mut sa = collect_segments(a, &everything);
    let mut sb = collect_segments(b, &everything);
    // Intersection (distance zero)?
    if ba.intersects(&bb) {
        let mut hit: Option<PointF> = None;
        any_pair(&mut sa, &mut sb, 0, |x, y| {
            if segments_intersect(x.0, x.1, y.0, y.1) {
                hit = Some(common_point(x.0, x.1, y.0, y.1));
                true
            } else {
                false
            }
        });
        if let Some(p) = hit {
            return Some(Closest {
                sq: SqDist::ZERO,
                a: p,
                b: p,
            });
        }
        if let Some(p) = containment(b, a).or_else(|| containment(a, b)) {
            return Some(Closest {
                sq: SqDist::ZERO,
                a: p.into(),
                b: p.into(),
            });
        }
    }
    // Nearest pair: branch-and-bound over a bounding-volume hierarchy of `b`, pruning with
    // the best distance found so far.
    let bvh = Bvh::build(&mut sb);
    let mut best = {
        let (s, p, q) = segment_segment(sa[0].0, sa[0].1, sb[0].0, sb[0].1);
        Closest { sq: s, a: p, b: q }
    };
    // Conservative integer radius bounding the current best distance.
    let radius =
        |c: &Closest| -> u128 { (libm::ceil(c.sq.distance_f64()) as u128).saturating_add(2) };
    let mut r = radius(&best);
    let mut stack: Vec<u32> = Vec::new();
    for x in &sa {
        stack.clear();
        stack.push(0);
        while let Some(ni) = stack.pop() {
            let node = &bvh.nodes[ni as usize];
            if rect_gap2(&x.2, &node.bbox) > r * r {
                continue;
            }
            if node.count > 0 {
                for y in &sb[node.first as usize..(node.first + node.count) as usize] {
                    if rect_gap2(&x.2, &y.2) > r * r {
                        continue;
                    }
                    let (s, p, q) = segment_segment(x.0, x.1, y.0, y.1);
                    if s < best.sq {
                        best = Closest { sq: s, a: p, b: q };
                        r = radius(&best);
                    }
                }
            } else {
                // Visit the nearer child first (pushed last).
                let (l, rr) = (node.first, node.first + 1);
                let gl = rect_gap2(&x.2, &bvh.nodes[l as usize].bbox);
                let gr = rect_gap2(&x.2, &bvh.nodes[rr as usize].bbox);
                if gl <= gr {
                    stack.push(rr);
                    stack.push(l);
                } else {
                    stack.push(l);
                    stack.push(rr);
                }
            }
        }
    }
    Some(best)
}

/// A bounding-volume hierarchy over segments (reordered in place): interior nodes have
/// `count == 0` and children at `first`, `first + 1`; leaves cover `first..first + count`.
struct Bvh {
    nodes: Vec<BvhNode>,
}

struct BvhNode {
    bbox: Rect,
    first: u32,
    count: u32,
}

impl Bvh {
    fn build(segs: &mut [(Point, Point, Rect)]) -> Bvh {
        let mut nodes = vec![BvhNode {
            bbox: bbox_of(segs),
            first: 0,
            count: segs.len() as u32,
        }];
        let mut stack = vec![0usize];
        while let Some(ni) = stack.pop() {
            let (first, count) = (nodes[ni].first as usize, nodes[ni].count as usize);
            if count <= 8 {
                continue;
            }
            let part = &mut segs[first..first + count];
            let b = nodes[ni].bbox;
            let mid = count / 2;
            // Split at the median centre along the longer side.
            if b.width() >= b.height() {
                part.select_nth_unstable_by_key(mid, |s| s.2.min.x as i128 + s.2.max.x as i128);
            } else {
                part.select_nth_unstable_by_key(mid, |s| s.2.min.y as i128 + s.2.max.y as i128);
            }
            let l = nodes.len();
            nodes.push(BvhNode {
                bbox: bbox_of(&part[..mid]),
                first: first as u32,
                count: mid as u32,
            });
            nodes.push(BvhNode {
                bbox: bbox_of(&part[mid..]),
                first: (first + mid) as u32,
                count: (count - mid) as u32,
            });
            nodes[ni] = BvhNode {
                bbox: b,
                first: l as u32,
                count: 0,
            };
            stack.push(l);
            stack.push(l + 1);
        }
        Bvh { nodes }
    }
}

fn bbox_of(segs: &[(Point, Point, Rect)]) -> Rect {
    segs.iter()
        .map(|s| s.2)
        .reduce(|a, b| a.union(&b))
        .unwrap_or(Rect {
            min: Point::new(0, 0),
            max: Point::new(0, 0),
        })
}

/// Exact squared minimum distance between two geometries (`None` when either is empty).
pub fn distance_sq<A: Geometry + ?Sized, B: Geometry + ?Sized>(a: &A, b: &B) -> Option<SqDist> {
    distance(a, b).map(|c| c.sq)
}

/// `true` when the distance between the two geometries is strictly less than `d`.
/// `false` for empty geometries and for coordinates outside
/// `±`[`MAX_COORD`](crate::MAX_COORD).
///
/// Exact, and much cheaper than [`distance`]: bounding boxes reject first, then only the
/// segments near the other geometry are compared, stopping at the first pair closer than
/// `d`; finally containment is tested. Always `false` for `d <= 0`.
///
/// ```
/// use polyclip::{distance_less_than, Path, Ring};
/// let pad = Ring::from([(0, 0), (100, 0), (100, 100), (0, 100)]);
/// let track = Path::from([(150, -50), (150, 150)]);
/// assert!(distance_less_than(&pad, &track, 51));
/// assert!(!distance_less_than(&pad, &track, 50));
/// ```
pub fn distance_less_than<A: Geometry + ?Sized, B: Geometry + ?Sized>(
    a: &A,
    b: &B,
    d: i64,
) -> bool {
    if d <= 0 {
        return false;
    }
    let (Some(ba), Some(bb)) = (a.bbox(), b.bbox()) else {
        return false;
    };
    if ![ba.min, ba.max, bb.min, bb.max]
        .iter()
        .all(|p| p.in_range())
    {
        return false;
    }
    let d2 = d as u128 * d as u128;
    if rect_gap2(&ba, &bb) >= d2 {
        return false;
    }
    let mut sa = collect_segments(a, &bb.expand(d));
    let mut sb = collect_segments(b, &ba.expand(d));
    if any_pair(&mut sa, &mut sb, d, |x, y| {
        rect_gap2(&x.2, &y.2) < d2 && segment_segment_lt(x.0, x.1, y.0, y.1, d2)
    }) {
        return true;
    }
    containment(b, a).is_some() || containment(a, b).is_some()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::{Path, Polygon, Ring};
    use crate::query::Segment;

    fn p(x: i64, y: i64) -> Point {
        Point::new(x, y)
    }
    fn sq(x0: i64, y0: i64, x1: i64, y1: i64) -> Ring {
        Ring::from([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])
    }

    #[test]
    fn basic() {
        let a = sq(0, 0, 10, 10);
        assert_eq!(
            distance(&a, &p(13, 14)).unwrap().sq.cmp_dist(5),
            Ordering::Equal
        );
        assert!(distance(&a, &p(5, 5)).unwrap().sq.is_zero());
        let h = Polygon::new(sq(0, 0, 10, 10), vec![sq(2, 2, 8, 8)]);
        // Point in the hole: distance 2 to the hole boundary.
        assert_eq!(
            distance(&h, &p(5, 4)).unwrap().sq.cmp_dist(2),
            Ordering::Equal
        );
        assert!(distance_less_than(&h, &p(5, 4), 3));
        assert!(!distance_less_than(&h, &p(5, 4), 2));
        // Diagonal: point (20, 20) to corner (10, 10): sqrt(200).
        let c = distance(&a, &p(20, 20)).unwrap();
        assert_eq!(c.sq, SqDist::int(200));
        assert_eq!(c.a, PointF::new(10.0, 10.0));
        // Non-integer interior projection: segment (0,0)-(10,3) to point (0,5): cross 50/sqrt(109).
        let s = Segment::new(p(0, 0), p(10, 3));
        let c = distance(&s, &p(0, 5)).unwrap();
        assert_eq!(c.sq, SqDist::ratio_sq(50, 109));
        assert!(distance_less_than(&s, &p(0, 5), 5));
        assert!(!distance_less_than(&s, &p(0, 5), 4));
        // Containment: small square inside big one.
        assert!(distance_less_than(&a, &sq(3, 3, 4, 4), 1));
        assert!(distance(&sq(3, 3, 4, 4), &a).unwrap().sq.is_zero());
        // Paths.
        let t = Path::from([(20, 0), (20, 10)]);
        assert_eq!(distance(&a, &t).unwrap().sq.cmp_dist(10), Ordering::Equal);
        assert!(!distance_less_than(&a, &t, 10));
        assert!(distance_less_than(&a, &t, 11));
    }

    #[test]
    fn huge_coordinates() {
        let m = crate::MAX_COORD;
        let s = Segment::new(p(-m, -m), p(m, m - 1));
        let q = p(-m, m);
        let c = distance(&s, &q).unwrap();
        let f = c.sq.distance_f64();
        assert!((f - (2.0f64 * m as f64) / 2f64.sqrt()).abs() < 1.0, "{f}");
        let di = f as i64;
        assert!(distance_less_than(&s, &q, di + 1));
        assert!(!distance_less_than(&s, &q, di));
    }

    #[test]
    fn brute_force_agreement() {
        let mut s: u64 = 99;
        let mut rnd = |m: i64| {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((s >> 33) as i64).rem_euclid(m)
        };
        for _ in 0..300 {
            let a: Path = (0..1 + rnd(6)).map(|_| p(rnd(100), rnd(100))).collect();
            let b: Path = (0..1 + rnd(6))
                .map(|_| p(rnd(100) + 50, rnd(100)))
                .collect();
            let c = distance(&a, &b).unwrap();
            // Brute force over all segment pairs.
            let segs = |x: &Path| -> Vec<(Point, Point)> {
                if x.len() == 1 {
                    vec![(x[0], x[0])]
                } else {
                    x.windows(2).map(|w| (w[0], w[1])).collect()
                }
            };
            let mut best: Option<SqDist> = None;
            for u in segs(&a) {
                for v in segs(&b) {
                    let s = segment_segment(u.0, u.1, v.0, v.1).0;
                    if best.is_none_or(|b| s < b) {
                        best = Some(s);
                    }
                }
            }
            assert_eq!(c.sq, best.unwrap());
            let df = c.sq.distance_f64();
            for d in [df.floor() as i64, df.ceil() as i64, df.ceil() as i64 + 1] {
                assert_eq!(
                    distance_less_than(&a, &b, d),
                    c.sq.cmp_dist(d.max(0) as u64) == Ordering::Less && d > 0
                );
            }
        }
    }
}
