//! Convex hull and Minkowski sum.

use crate::boolean::{FillRule, union_all};
use crate::error::{Error, Result, check_point};
use crate::geom::{Point, Polygon, PolygonSet, Ring};
use crate::predicates::orient;
use crate::query::ring_area2;

/// Convex hull of a point set (Andrew's monotone chain, exact).
///
/// Returns a canonical counter-clockwise ring starting at the lexicographically smallest
/// point, without collinear vertices. Degenerate inputs give degenerate rings: one point for
/// a single (repeated) point, the two extreme points when all points are collinear, and an
/// empty ring for no points. Coordinates outside `±`[`MAX_COORD`](crate::MAX_COORD) are an
/// error.
///
/// ```
/// use polyclip::{convex_hull, Point, Ring};
/// let pts = [(0, 0), (4, 0), (2, 1), (4, 4), (0, 4), (2, 2)].map(Point::from);
/// assert_eq!(convex_hull(pts).unwrap(), Ring::from([(0, 0), (4, 0), (4, 4), (0, 4)]));
/// ```
pub fn convex_hull(points: impl IntoIterator<Item = Point>) -> Result<Ring> {
    let mut p: Vec<Point> = points.into_iter().collect();
    for &q in &p {
        check_point(q)?;
    }
    p.sort_unstable();
    p.dedup();
    if p.len() <= 2 {
        return Ok(Ring(p));
    }
    let mut h: Vec<Point> = Vec::with_capacity(p.len() + 1);
    // Lower hull.
    for &q in &p {
        while h.len() >= 2 && orient(h[h.len() - 2], h[h.len() - 1], q) <= 0 {
            h.pop();
        }
        h.push(q);
    }
    // Upper hull.
    let lower = h.len() + 1;
    for &q in p.iter().rev().skip(1) {
        while h.len() >= lower && orient(h[h.len() - 2], h[h.len() - 1], q) <= 0 {
            h.pop();
        }
        h.push(q);
    }
    h.pop();
    Ok(Ring(h))
}

/// Convex hull of every vertex of a geometry (rings, polygons, paths, ...). See
/// [`convex_hull`].
pub fn convex_hull_of<G: crate::query::Geometry + ?Sized>(g: &G) -> Result<Ring> {
    let mut pts = Vec::new();
    g.visit_segments(&mut |a, b| {
        pts.push(a);
        pts.push(b);
    });
    convex_hull(pts)
}

/// Convex and simple: all turns have the sign of the area, and the ring winds exactly once
/// (its area equals its hull's, which rules out stars whose turns all agree).
fn is_convex(r: &[Point]) -> bool {
    let n = r.len();
    if n < 3 {
        return false;
    }
    let a = ring_area2(r);
    let s = a.signum();
    if s == 0 || !(0..n).all(|i| orient(r[i], r[(i + 1) % n], r[(i + 2) % n]).signum() * s >= 0) {
        return false;
    }
    convex_hull(r.iter().copied()).is_ok_and(|h| h.signed_area2() == a.abs())
}

fn add(a: Point, b: Point) -> Result<Point> {
    let p = Point::new(a.x + b.x, a.y + b.y);
    check_point(p)?;
    Ok(p)
}

/// Minkowski sum `A ⊕ B = { a + b : a ∈ A, b ∈ B }` of two polygons (with holes), exact.
///
/// Computed as the union of the parallelograms `e ⊕ f` for every edge pair, the copies of
/// `B` translated to every vertex of `A`, and `A` translated by one vertex of `B` — an
/// identity that holds for any connected `B`. When both polygons are convex and hole-free,
/// the convex hull of pairwise vertex sums is returned directly. Inputs are normalized by
/// orientation; the result is canonical. Sums outside `±MAX_COORD` are an error.
///
/// ```
/// use polyclip::{minkowski_sum, Polygon, Ring};
/// let a = Polygon::from(Ring::from([(0, 0), (10, 0), (10, 10), (0, 10)]));
/// let b = Polygon::from(Ring::from([(-1, -1), (1, -1), (1, 1), (-1, 1)]));
/// let s = minkowski_sum(&a, &b).unwrap();
/// assert_eq!(s[0].outer, Ring::from([(-1, -1), (11, -1), (11, 11), (-1, 11)]));
/// ```
pub fn minkowski_sum(a: &Polygon, b: &Polygon) -> Result<PolygonSet> {
    for p in a.rings().chain(b.rings()).flat_map(|r| r.iter()) {
        check_point(*p)?;
    }
    if a.outer.len() < 3 || b.outer.len() < 3 {
        return Err(Error::InvalidParameter(
            "minkowski_sum needs polygons with at least 3 vertices",
        ));
    }
    if a.holes.is_empty() && b.holes.is_empty() && is_convex(&a.outer) && is_convex(&b.outer) {
        let mut sums = Vec::with_capacity(a.outer.len() * b.outer.len());
        for &p in a.outer.iter() {
            for &q in b.outer.iter() {
                sums.push(add(p, q)?);
            }
        }
        let h = convex_hull(sums)?;
        return union_all(&h, FillRule::NonZero);
    }
    // Normalize each operand to a clean region first.
    let na = union_all(a, FillRule::NonZero)?;
    let nb = union_all(b, FillRule::NonZero)?;
    let mut out: PolygonSet = Vec::new();
    for pb in &nb {
        let mut pieces: Vec<Ring> = Vec::new();
        let b0 = pb.outer[0];
        for pa in &na {
            // A + b0 (with its holes).
            for r in pa.rings() {
                pieces.push(r.iter().map(|&p| add(p, b0)).collect::<Result<Ring>>()?);
            }
            for ra in pa.rings() {
                let n = ra.len();
                for i in 0..n {
                    let (p0, p1) = (ra[i], ra[(i + 1) % n]);
                    // B + p0 (with its holes).
                    for rb in pb.rings() {
                        pieces.push(rb.iter().map(|&q| add(q, p0)).collect::<Result<Ring>>()?);
                    }
                    // Edge-pair parallelograms.
                    for rb in pb.rings() {
                        let m = rb.len();
                        for j in 0..m {
                            let (q0, q1) = (rb[j], rb[(j + 1) % m]);
                            let quad = [add(p0, q0)?, add(p1, q0)?, add(p1, q1)?, add(p0, q1)?];
                            let mut quad = Ring(quad.to_vec());
                            if quad.signed_area2() < 0 {
                                quad.reverse_orientation();
                            }
                            pieces.push(quad);
                        }
                    }
                }
            }
        }
        // Holes of translated copies are clockwise: NonZero union of the pieces is the sum.
        out.extend(union_all(&pieces, FillRule::Positive)?);
    }
    union_all(&out, FillRule::NonZero)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{FillRule, Op, boolean};

    fn sq(x0: i64, y0: i64, x1: i64, y1: i64) -> Ring {
        Ring::from([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])
    }

    #[test]
    fn hull_degenerate() {
        assert!(convex_hull(Vec::new()).unwrap().is_empty());
        assert_eq!(convex_hull([Point::new(1, 1); 3]).unwrap().len(), 1);
        assert!(
            convex_hull([
                Point::new(i64::MIN, 0),
                Point::new(i64::MAX, 1),
                Point::new(0, 5)
            ])
            .is_err()
        );
        assert_eq!(
            convex_hull([(0, 0), (1, 1), (2, 2)].map(Point::from)).unwrap(),
            Ring::from([(0, 0), (2, 2)])
        );
    }

    #[test]
    fn minkowski_nonconvex() {
        // L shape plus a small square: equals the L offset with square joins by 1.
        let l = Polygon::from(Ring::from([
            (0, 0),
            (20, 0),
            (20, 10),
            (10, 10),
            (10, 20),
            (0, 20),
        ]));
        let b = Polygon::from(sq(-1, -1, 1, 1));
        let s = minkowski_sum(&l, &b).unwrap();
        let e = crate::offset(
            &l,
            1,
            crate::Join::Miter { limit: 2.0 },
            crate::ArcTol::new(1, crate::Side::Inside),
        )
        .unwrap();
        assert_eq!(s, e);
        // Polygon with a hole that closes up.
        let h = Polygon::new(
            sq(0, 0, 10, 10),
            vec![Ring::from([(4, 4), (4, 6), (6, 6), (6, 4)])],
        );
        let s = minkowski_sum(&h, &b).unwrap();
        assert_eq!(s, vec![Polygon::from(sq(-1, -1, 11, 11))]);
        let small = Polygon::from(Ring::from([(0, 0), (1, 0), (1, 1), (0, 1)]));
        let s = minkowski_sum(&h, &small).unwrap();
        assert_eq!(s[0].holes.len(), 1);
        assert_eq!(s[0].holes[0], Ring::from([(5, 5), (5, 6), (6, 6), (6, 5)]));
        // Commutative.
        let s1 = minkowski_sum(&l, &h).unwrap();
        let s2 = minkowski_sum(&h, &l).unwrap();
        assert!(
            boolean(Op::Xor, &s1, &s2, FillRule::NonZero)
                .unwrap()
                .is_empty()
        );
    }
}
