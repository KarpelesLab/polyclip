//! Measurements and exact point/region queries.
//!
//! All predicates here are exact. Regions are closed sets: a point on a boundary is in the
//! region, and two regions touching at a single point intersect.

use crate::arrangement::{Arrangement, InEdge, Noding};
use crate::geom::{Path, Point, PointF, PolyTree, Polygon, Rect, Ring};
use crate::predicates::{on_segment, orient, segments_intersect};

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
        s += orient(o, pts[i], pts[i + 1]);
    }
    s
}

/// Where a point lies relative to a geometry.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Location {
    /// In the interior.
    Inside,
    /// Outside.
    Outside,
    /// On the boundary (for linear geometries: on the geometry itself).
    OnBoundary,
}

/// A line segment between two integer points.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Segment {
    /// Start point.
    pub a: Point,
    /// End point.
    pub b: Point,
}

impl Segment {
    /// Creates a segment.
    pub const fn new(a: Point, b: Point) -> Self {
        Segment { a, b }
    }
}

/// Winding number of `pts` (a closed ring) around `p`, or `None` when `p` is on the ring.
pub fn ring_winding(pts: &[Point], p: Point) -> Option<i32> {
    let n = pts.len();
    let mut w = 0i32;
    for i in 0..n {
        let a = pts[i];
        let b = pts[if i + 1 == n { 0 } else { i + 1 }];
        if a.y <= p.y {
            if b.y > p.y {
                let o = orient(a, b, p);
                if o > 0 {
                    w += 1;
                } else if o == 0 {
                    return None;
                }
            } else if b.y == p.y && on_segment(a, b, p) {
                return None;
            }
        } else if b.y <= p.y {
            let o = orient(a, b, p);
            if o < 0 {
                w -= 1;
            } else if o == 0 {
                return None;
            }
        }
    }
    Some(w)
}

/// Location of `p` relative to the region bounded by a ring (non-zero winding rule; for
/// simple rings this is the usual interior).
pub fn locate_in_ring(pts: &[Point], p: Point) -> Location {
    match ring_winding(pts, p) {
        None => Location::OnBoundary,
        Some(0) => Location::Outside,
        Some(_) => Location::Inside,
    }
}

/// A geometry that queries ([`locate`], [`intersects`], [`contains`], distances) accept.
///
/// Areal geometries ([`Ring`], [`Polygon`], polygon sets, [`PolyTree`]) are closed regions;
/// linear ones ([`Segment`], [`Path`]) and [`Point`] are just their point sets.
pub trait Geometry {
    /// Bounding box, `None` when empty.
    fn bbox(&self) -> Option<Rect>;
    /// `true` for regions (rings, polygons), `false` for points and linear geometries.
    fn is_areal(&self) -> bool;
    /// Calls `f` for every boundary segment (areal) or every segment (linear). A point
    /// geometry yields one degenerate segment `(p, p)`; so does a single-vertex path.
    fn visit_segments(&self, f: &mut dyn FnMut(Point, Point));
    /// Location of `p` relative to the geometry.
    fn locate(&self, p: Point) -> Location;
    /// Some vertex of the geometry, if not empty.
    fn any_point(&self) -> Option<Point>;
}

impl Geometry for Point {
    fn bbox(&self) -> Option<Rect> {
        Some(Rect {
            min: *self,
            max: *self,
        })
    }
    fn is_areal(&self) -> bool {
        false
    }
    fn visit_segments(&self, f: &mut dyn FnMut(Point, Point)) {
        f(*self, *self)
    }
    fn locate(&self, p: Point) -> Location {
        if p == *self {
            Location::OnBoundary
        } else {
            Location::Outside
        }
    }
    fn any_point(&self) -> Option<Point> {
        Some(*self)
    }
}

impl Geometry for Segment {
    fn bbox(&self) -> Option<Rect> {
        Some(Rect::new(self.a, self.b))
    }
    fn is_areal(&self) -> bool {
        false
    }
    fn visit_segments(&self, f: &mut dyn FnMut(Point, Point)) {
        f(self.a, self.b)
    }
    fn locate(&self, p: Point) -> Location {
        if on_segment(self.a, self.b, p) {
            Location::OnBoundary
        } else {
            Location::Outside
        }
    }
    fn any_point(&self) -> Option<Point> {
        Some(self.a)
    }
}

impl Geometry for Path {
    fn bbox(&self) -> Option<Rect> {
        Rect::of_points(self.0.iter())
    }
    fn is_areal(&self) -> bool {
        false
    }
    fn visit_segments(&self, f: &mut dyn FnMut(Point, Point)) {
        if self.0.len() == 1 {
            f(self.0[0], self.0[0]);
        }
        for w in self.0.windows(2) {
            f(w[0], w[1]);
        }
    }
    fn locate(&self, p: Point) -> Location {
        let on = if self.0.len() == 1 {
            self.0[0] == p
        } else {
            self.0.windows(2).any(|w| on_segment(w[0], w[1], p))
        };
        if on {
            Location::OnBoundary
        } else {
            Location::Outside
        }
    }
    fn any_point(&self) -> Option<Point> {
        self.0.first().copied()
    }
}

fn visit_ring(pts: &[Point], f: &mut dyn FnMut(Point, Point)) {
    let n = pts.len();
    for i in 0..n {
        f(pts[i], pts[(i + 1) % n]);
    }
}

impl Geometry for Ring {
    fn bbox(&self) -> Option<Rect> {
        Rect::of_points(self.0.iter())
    }
    fn is_areal(&self) -> bool {
        true
    }
    fn visit_segments(&self, f: &mut dyn FnMut(Point, Point)) {
        visit_ring(&self.0, f)
    }
    fn locate(&self, p: Point) -> Location {
        locate_in_ring(&self.0, p)
    }
    fn any_point(&self) -> Option<Point> {
        self.0.first().copied()
    }
}

/// Location of `p` relative to a polygon (inside its outer ring and outside all holes).
pub fn locate_in_polygon(poly: &Polygon, p: Point) -> Location {
    if let Some(b) = poly.outer.bbox() {
        if !b.contains_point(p) {
            return Location::Outside;
        }
    }
    match locate_in_ring(&poly.outer.0, p) {
        Location::Inside => {}
        other => return other,
    }
    for h in &poly.holes {
        match locate_in_ring(&h.0, p) {
            Location::Outside => {}
            Location::Inside => return Location::Outside,
            Location::OnBoundary => return Location::OnBoundary,
        }
    }
    Location::Inside
}

impl Geometry for Polygon {
    fn bbox(&self) -> Option<Rect> {
        self.outer.bbox()
    }
    fn is_areal(&self) -> bool {
        true
    }
    fn visit_segments(&self, f: &mut dyn FnMut(Point, Point)) {
        for r in self.rings() {
            visit_ring(&r.0, f)
        }
    }
    fn locate(&self, p: Point) -> Location {
        locate_in_polygon(self, p)
    }
    fn any_point(&self) -> Option<Point> {
        self.outer.0.first().copied()
    }
}

impl Geometry for [Polygon] {
    fn bbox(&self) -> Option<Rect> {
        self.iter()
            .filter_map(|p| p.bbox())
            .reduce(|a, b| a.union(&b))
    }
    fn is_areal(&self) -> bool {
        true
    }
    fn visit_segments(&self, f: &mut dyn FnMut(Point, Point)) {
        for p in self {
            p.visit_segments(f)
        }
    }
    fn locate(&self, p: Point) -> Location {
        let mut res = Location::Outside;
        for poly in self {
            match poly.locate(p) {
                Location::Inside => return Location::Inside,
                Location::OnBoundary => res = Location::OnBoundary,
                Location::Outside => {}
            }
        }
        res
    }
    fn any_point(&self) -> Option<Point> {
        self.iter().find_map(|p| p.any_point())
    }
}

impl Geometry for Vec<Polygon> {
    fn bbox(&self) -> Option<Rect> {
        self.as_slice().bbox()
    }
    fn is_areal(&self) -> bool {
        true
    }
    fn visit_segments(&self, f: &mut dyn FnMut(Point, Point)) {
        self.as_slice().visit_segments(f)
    }
    fn locate(&self, p: Point) -> Location {
        self.as_slice().locate(p)
    }
    fn any_point(&self) -> Option<Point> {
        self.as_slice().any_point()
    }
}

impl Geometry for PolyTree {
    fn bbox(&self) -> Option<Rect> {
        self.nodes
            .iter()
            .filter_map(|n| n.ring.bbox())
            .reduce(|a, b| a.union(&b))
    }
    fn is_areal(&self) -> bool {
        true
    }
    fn visit_segments(&self, f: &mut dyn FnMut(Point, Point)) {
        for n in &self.nodes {
            visit_ring(&n.ring.0, f)
        }
    }
    fn locate(&self, p: Point) -> Location {
        // Count enclosing rings: odd means inside (valid nesting alternates outer/hole).
        let mut count = 0usize;
        for n in &self.nodes {
            match locate_in_ring(&n.ring.0, p) {
                Location::OnBoundary => return Location::OnBoundary,
                Location::Inside => count += 1,
                Location::Outside => {}
            }
        }
        if count % 2 == 1 {
            Location::Inside
        } else {
            Location::Outside
        }
    }
    fn any_point(&self) -> Option<Point> {
        self.nodes.first().and_then(|n| n.ring.0.first().copied())
    }
}

impl<T: Geometry + ?Sized> Geometry for &T {
    fn bbox(&self) -> Option<Rect> {
        (**self).bbox()
    }
    fn is_areal(&self) -> bool {
        (**self).is_areal()
    }
    fn visit_segments(&self, f: &mut dyn FnMut(Point, Point)) {
        (**self).visit_segments(f)
    }
    fn locate(&self, p: Point) -> Location {
        (**self).locate(p)
    }
    fn any_point(&self) -> Option<Point> {
        (**self).any_point()
    }
}

/// Location of `p` relative to `g`.
pub fn locate<G: Geometry + ?Sized>(g: &G, p: Point) -> Location {
    g.locate(p)
}

/// Segments of `g` whose bounding box meets `filter`, with their bounding boxes.
pub(crate) fn collect_segments<G: Geometry + ?Sized>(
    g: &G,
    filter: &Rect,
) -> Vec<(Point, Point, Rect)> {
    let mut v = Vec::new();
    g.visit_segments(&mut |a, b| {
        let r = Rect::new(a, b);
        if r.intersects(filter) {
            v.push((a, b, r));
        }
    });
    v
}

/// Calls `f` on every pair `(sa[i], sb[j])` whose bounding boxes, grown by `margin`,
/// overlap, stopping as soon as `f` returns `true`. Returns whether it stopped.
pub(crate) fn any_pair(
    sa: &mut [(Point, Point, Rect)],
    sb: &mut [(Point, Point, Rect)],
    margin: i64,
    mut f: impl FnMut(&(Point, Point, Rect), &(Point, Point, Rect)) -> bool,
) -> bool {
    if sa.len() * sb.len() <= 256 {
        for x in sa.iter() {
            let rx = x.2.expand(margin);
            for y in sb.iter() {
                if rx.intersects(&y.2) && f(x, y) {
                    return true;
                }
            }
        }
        return false;
    }
    // Sweep-and-prune along x: merge both lists by min x, keeping for each side the items
    // whose x-extent still reaches the sweep position.
    sa.sort_unstable_by_key(|s| s.2.min.x);
    sb.sort_unstable_by_key(|s| s.2.min.x);
    let (mut i, mut j) = (0usize, 0usize);
    let mut act_a: Vec<usize> = Vec::new();
    let mut act_b: Vec<usize> = Vec::new();
    while i < sa.len() || j < sb.len() {
        let take_a = j >= sb.len() || (i < sa.len() && sa[i].2.min.x <= sb[j].2.min.x);
        if take_a {
            let x = &sa[i];
            let rx = x.2.expand(margin);
            act_b.retain(|&k| sb[k].2.max.x >= rx.min.x);
            for &k in &act_b {
                if rx.intersects(&sb[k].2) && f(x, &sb[k]) {
                    return true;
                }
            }
            act_a.push(i);
            i += 1;
        } else {
            let y = &sb[j];
            let ry = y.2.expand(margin);
            act_a.retain(|&k| sa[k].2.max.x >= ry.min.x);
            for &k in &act_a {
                if ry.intersects(&sa[k].2) && f(&sa[k], y) {
                    return true;
                }
            }
            act_b.push(j);
            j += 1;
        }
    }
    false
}

/// `true` when the two geometries share at least one point (closed sets).
///
/// Exact. Bounding boxes reject early; then boundaries are tested for any contact, and
/// finally one geometry is tested for lying inside the other.
pub fn intersects<A: Geometry + ?Sized, B: Geometry + ?Sized>(a: &A, b: &B) -> bool {
    let (Some(ba), Some(bb)) = (a.bbox(), b.bbox()) else {
        return false;
    };
    if !ba.intersects(&bb) {
        return false;
    }
    let mut sa = collect_segments(a, &bb);
    let mut sb = collect_segments(b, &ba);
    if any_pair(&mut sa, &mut sb, 0, |x, y| {
        segments_intersect(x.0, x.1, y.0, y.1)
    }) {
        return true;
    }
    if b.is_areal() {
        if let Some(p) = a.any_point() {
            if b.locate(p) != Location::Outside {
                return true;
            }
        }
    }
    if a.is_areal() {
        if let Some(p) = b.any_point() {
            if a.locate(p) != Location::Outside {
                return true;
            }
        }
    }
    false
}

/// `true` when every point of `b` belongs to `a` (closed sets), exactly.
///
/// `a` must be areal for this to be meaningful; for a linear or point `a` this tests whether
/// `b` is a point lying on `a`. An empty `b` is contained in anything. Inputs are expected
/// to be valid (non-self-crossing); for self-crossing input the result is unspecified but
/// the function never panics.
pub fn contains<A: Geometry + ?Sized, B: Geometry + ?Sized>(a: &A, b: &B) -> bool {
    let Some(bb) = b.bbox() else {
        return true;
    };
    let Some(ba) = a.bbox() else {
        return false;
    };
    if !ba.contains_rect(&bb) {
        return false;
    }
    if !a.is_areal() {
        let mut pts = Vec::new();
        let mut linear = false;
        b.visit_segments(&mut |p, q| {
            if p == q { pts.push(p) } else { linear = true }
        });
        return !linear && !b.is_areal() && pts.iter().all(|&p| a.locate(p) != Location::Outside);
    }
    let mut edges: Vec<InEdge> = Vec::new();
    a.visit_segments(&mut |p, q| {
        if p != q {
            edges.push(InEdge {
                a: p,
                b: q,
                tag: 0,
                operand: 0,
            });
        }
    });
    let b_areal = b.is_areal();
    let mut isolated = Vec::new();
    b.visit_segments(&mut |p, q| {
        if p == q {
            isolated.push(p);
        } else {
            edges.push(InEdge {
                a: p,
                b: q,
                tag: 0,
                operand: if b_areal { 1 } else { 2 },
            });
        }
    });
    if isolated.iter().any(|&p| a.locate(p) == Location::Outside) {
        return false;
    }
    let Ok(arr) = Arrangement::build(&edges, Noding::Exact) else {
        return false;
    };
    let in_a = |w: [i32; 2]| w[0] & 1 != 0;
    let in_b = |w: [i32; 2]| w[1] & 1 != 0;
    for k in 0..arr.edges.len() {
        let (wb, wa) = arr.sides(k);
        if arr.edges[k].open.is_some() {
            if !in_a(wb) && !in_a(wa) {
                return false;
            }
        } else if (in_b(wb) && !in_a(wb)) || (in_b(wa) && !in_a(wa)) {
            return false;
        }
    }
    true
}

/// Twice the signed area of a geometry's region (exact): the sum over its rings.
pub fn area2<G: Geometry + ?Sized>(g: &G) -> i128 {
    if !g.is_areal() {
        return 0;
    }
    // Rings are visited edge by edge; the shoelace sum is translation invariant per closed
    // ring, so summing cross products about the origin works ring by ring.
    let mut s: i128 = 0;
    g.visit_segments(&mut |a, b| s += a.x as i128 * b.y as i128 - b.x as i128 * a.y as i128);
    s
}

/// Area centroid of a region given as rings (holes clockwise subtract), or `None` when the
/// area is zero. Computed in `f64` from exact per-edge terms, relative to the first vertex.
pub fn centroid<G: Geometry + ?Sized>(g: &G) -> Option<PointF> {
    if !g.is_areal() {
        return None;
    }
    let o = g.any_point()?;
    let mut a2: i128 = 0;
    let (mut cx, mut cy) = (0f64, 0f64);
    g.visit_segments(&mut |p, q| {
        let px = (p.x - o.x) as i128;
        let py = (p.y - o.y) as i128;
        let qx = (q.x - o.x) as i128;
        let qy = (q.y - o.y) as i128;
        let c = px * qy - qx * py;
        a2 += c;
        cx += (px + qx) as f64 * c as f64;
        cy += (py + qy) as f64 * c as f64;
    });
    if a2 == 0 {
        return None;
    }
    let d = 3.0 * a2 as f64;
    Some(PointF::new(o.x as f64 + cx / d, o.y as f64 + cy / d))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sq(x0: i64, y0: i64, x1: i64, y1: i64) -> Ring {
        Ring::from([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])
    }
    fn p(x: i64, y: i64) -> Point {
        Point::new(x, y)
    }

    #[test]
    fn locate_ring() {
        let r = sq(0, 0, 10, 10);
        assert_eq!(locate(&r, p(5, 5)), Location::Inside);
        assert_eq!(locate(&r, p(10, 5)), Location::OnBoundary);
        assert_eq!(locate(&r, p(0, 0)), Location::OnBoundary);
        assert_eq!(locate(&r, p(5, 10)), Location::OnBoundary);
        assert_eq!(locate(&r, p(11, 5)), Location::Outside);
        assert_eq!(locate(&r, p(-1, 10)), Location::Outside);
        let poly = Polygon::new(sq(0, 0, 10, 10), vec![sq(2, 2, 4, 4)]);
        assert_eq!(locate(&poly, p(3, 3)), Location::Outside);
        assert_eq!(locate(&poly, p(2, 3)), Location::OnBoundary);
        assert_eq!(locate(&poly, p(5, 5)), Location::Inside);
    }

    #[test]
    fn intersects_contains() {
        let a = sq(0, 0, 10, 10);
        let b = sq(10, 10, 20, 20);
        let c = sq(2, 2, 4, 4);
        assert!(intersects(&a, &b));
        assert!(!intersects(&a, &sq(11, 0, 20, 20)));
        assert!(intersects(&a, &c));
        assert!(contains(&a, &c));
        assert!(!contains(&c, &a));
        assert!(contains(&a, &a));
        assert!(!contains(&a, &b));
        let holed = Polygon::new(sq(0, 0, 10, 10), vec![sq(2, 2, 4, 4)]);
        assert!(!contains(&holed, &c));
        assert!(!contains(&holed, &sq(1, 1, 5, 5)));
        assert!(contains(&holed, &sq(4, 4, 6, 6)));
        assert!(!contains(&holed, &Segment::new(p(0, 0), p(10, 10))));
        assert!(contains(&holed, &Segment::new(p(0, 0), p(2, 2))));
        assert!(contains(&holed, &Segment::new(p(0, 0), p(10, 0))));
        assert!(contains(&a, &p(10, 10)));
        assert!(!contains(&a, &p(11, 10)));
        assert!(!intersects(&holed, &p(3, 3)));
        assert!(intersects(&holed, &p(2, 3)));
        // Square in hole exactly equal to the hole: not contained, touching boundary.
        assert!(!contains(&holed, &sq(2, 2, 4, 4)));
        assert!(intersects(&holed, &sq(2, 2, 4, 4)));
    }

    #[test]
    fn areas() {
        assert_eq!(area2(&sq(0, 0, 10, 10)), 200);
        let holed = Polygon::new(
            sq(0, 0, 10, 10),
            vec![Ring::from([(2, 2), (2, 4), (4, 4), (4, 2)])],
        );
        assert_eq!(area2(&holed), 2 * 96);
        let c = centroid(&sq(0, 0, 10, 4)).unwrap();
        assert!((c.x - 5.0).abs() < 1e-12 && (c.y - 2.0).abs() < 1e-12);
    }
}
