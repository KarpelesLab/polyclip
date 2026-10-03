//! Validity and canonical-form checks.

use crate::arrangement::{Arrangement, InEdge};
use crate::geom::{Point, PointF, Polygon};
use crate::node::node_exact;
use crate::predicates::{crossing_f64, orient};
use crate::query::{Location, locate_in_ring, ring_area2};
use core::fmt;

/// Identifies a ring: `polygon` index in the set (0 for a single polygon) and `ring` index
/// within the polygon (0 = outer, `1 + i` = hole `i`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct RingId {
    /// Polygon index.
    pub polygon: usize,
    /// Ring index (0 = outer ring).
    pub ring: usize,
}

/// Why a polygon (set) is invalid.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum ValidityError {
    /// A coordinate is outside `±MAX_COORD`.
    CoordinateOutOfRange(Point),
    /// A ring has fewer than three vertices.
    TooFewVertices(RingId),
    /// Two consecutive vertices are equal (a zero-length edge).
    DuplicateVertex(RingId, Point),
    /// A ring has zero area.
    ZeroArea(RingId),
    /// Two edges cross at a point interior to both.
    SelfIntersection(PointF),
    /// A ring touches itself (passes twice through a point).
    SelfTouch(RingId, Point),
    /// Two edges overlap along a segment.
    OverlappingEdges(Point, Point),
    /// A hole is not inside its outer ring, holes overlap, or polygons overlap. The point is
    /// an endpoint of an edge bordering the offending region.
    InvalidNesting(Point),
    /// The rings touch in a way that splits the polygon's interior into several pieces.
    DisconnectedInterior(Point),
    /// Orientation does not match the canonical convention (outer CCW, holes CW).
    WrongOrientation(RingId),
    /// Not in canonical form (start vertex, collinear vertex or ordering).
    NotCanonical(&'static str),
}

impl fmt::Display for ValidityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use ValidityError::*;
        match self {
            CoordinateOutOfRange(p) => write!(f, "coordinate out of range at ({}, {})", p.x, p.y),
            TooFewVertices(r) => {
                write!(f, "ring {}/{} has fewer than 3 vertices", r.polygon, r.ring)
            }
            DuplicateVertex(r, p) => write!(
                f,
                "ring {}/{} repeats vertex ({}, {})",
                r.polygon, r.ring, p.x, p.y
            ),
            ZeroArea(r) => write!(f, "ring {}/{} has zero area", r.polygon, r.ring),
            SelfIntersection(p) => write!(f, "self-intersection at ({}, {})", p.x, p.y),
            SelfTouch(r, p) => write!(
                f,
                "ring {}/{} touches itself at ({}, {})",
                r.polygon, r.ring, p.x, p.y
            ),
            OverlappingEdges(a, b) => write!(
                f,
                "overlapping edges along ({}, {})-({}, {})",
                a.x, a.y, b.x, b.y
            ),
            InvalidNesting(p) => write!(f, "invalid nesting near ({}, {})", p.x, p.y),
            DisconnectedInterior(p) => write!(f, "interior disconnected at ({}, {})", p.x, p.y),
            WrongOrientation(r) => write!(
                f,
                "ring {}/{} has non-canonical orientation",
                r.polygon, r.ring
            ),
            NotCanonical(s) => write!(f, "not canonical: {s}"),
        }
    }
}

impl std::error::Error for ValidityError {}

/// Checks that a polygon is valid: every ring simple with at least three vertices and
/// non-zero area, holes strictly inside the outer ring and with disjoint interiors, rings
/// touching only at isolated points without disconnecting the interior. Orientation is not
/// checked (see [`check_canonical`]).
pub fn validate(poly: &Polygon) -> Result<(), ValidityError> {
    validate_set(core::slice::from_ref(poly))
}

/// Checks a polygon set: every polygon valid (see [`validate`]) and polygon interiors
/// pairwise disjoint (touching at isolated points is allowed, sharing edges is not).
pub fn validate_set(polys: &[Polygon]) -> Result<(), ValidityError> {
    // Per-ring checks, and edges with each ring oriented by its role.
    let mut edges: Vec<InEdge> = Vec::new();
    let mut ring_of_edge: Vec<(usize, usize)> = Vec::new();
    for (pi, poly) in polys.iter().enumerate() {
        for (ri, ring) in poly.rings().enumerate() {
            let id = RingId {
                polygon: pi,
                ring: ri,
            };
            let pts = &ring.0;
            if let Some(&p) = pts.iter().find(|p| !p.in_range()) {
                return Err(ValidityError::CoordinateOutOfRange(p));
            }
            if pts.len() < 3 {
                return Err(ValidityError::TooFewVertices(id));
            }
            for i in 0..pts.len() {
                if pts[i] == pts[(i + 1) % pts.len()] {
                    return Err(ValidityError::DuplicateVertex(id, pts[i]));
                }
            }
            let a = ring_area2(pts);
            if a == 0 {
                return Err(ValidityError::ZeroArea(id));
            }
            // Orient outer rings CCW and holes CW for the winding check.
            let flip = (ri == 0) != (a > 0);
            let n = pts.len();
            for i in 0..n {
                let (p, q) = (pts[i], pts[(i + 1) % n]);
                let (p, q) = if flip { (q, p) } else { (p, q) };
                edges.push(InEdge {
                    a: p,
                    b: q,
                    tag: ((pi as u64) << 32) | ri as u64,
                    operand: 0,
                });
                ring_of_edge.push((pi, ri));
            }
        }
    }
    let segs: Vec<(Point, Point)> = edges.iter().map(|e| (e.a, e.b)).collect();
    let frags = match node_exact(&segs) {
        Ok(f) => f,
        Err(c) => {
            let (a, b) = (segs[c.i as usize], segs[c.j as usize]);
            let (x, y) = crossing_f64(a.0, a.1, b.0, b.1);
            return Err(ValidityError::SelfIntersection(PointF::new(x, y)));
        }
    };
    drop(segs);
    // Overlapping edges.
    let mut fr: Vec<(Point, Point)> = frags
        .iter()
        .map(|f| if f.a < f.b { (f.a, f.b) } else { (f.b, f.a) })
        .collect();
    fr.sort_unstable();
    for w in fr.windows(2) {
        if w[0] == w[1] {
            return Err(ValidityError::OverlappingEdges(w[0].0, w[0].1));
        }
    }
    drop(fr);
    // Self-touch (a ring through a point more than once) and the touch graph.
    let mut vr: Vec<(Point, (usize, usize))> = Vec::with_capacity(frags.len() * 2);
    for f in &frags {
        let r = ring_of_edge[f.src as usize];
        vr.push((f.a, r));
        vr.push((f.b, r));
    }
    vr.sort_unstable();
    // Union-find over rings and touch points (per polygon): a cycle means the touching
    // rings cut the interior into pieces.
    let mut uf = UnionFind::new(0);
    let mut ring_node: std::collections::HashMap<(usize, usize), usize> =
        std::collections::HashMap::new();
    let mut rings_here: Vec<(usize, usize)> = Vec::new();
    let mut i = 0;
    while i < vr.len() {
        let p = vr[i].0;
        let mut j = i;
        rings_here.clear();
        while j < vr.len() && vr[j].0 == p {
            let r = vr[j].1;
            let mut k = j;
            while k < vr.len() && vr[k].0 == p && vr[k].1 == r {
                k += 1;
            }
            // A simple ring passes through each of its vertices exactly once (two edges).
            if k - j > 2 {
                return Err(ValidityError::SelfTouch(
                    RingId {
                        polygon: r.0,
                        ring: r.1,
                    },
                    p,
                ));
            }
            rings_here.push(r);
            j = k;
        }
        let mut a = 0;
        while a < rings_here.len() {
            let mut b = a;
            while b < rings_here.len() && rings_here[b].0 == rings_here[a].0 {
                b += 1;
            }
            if b - a >= 2 {
                let pn = uf.add();
                for r in &rings_here[a..b] {
                    let rn = *ring_node.entry(*r).or_insert_with(|| uf.add());
                    if !uf.union(pn, rn) {
                        return Err(ValidityError::DisconnectedInterior(p));
                    }
                }
            }
            a = b;
        }
        i = j;
    }
    drop(vr);
    let arr = Arrangement::from_frags(&edges, frags);
    // Winding numbers must be 0 or 1 everywhere.
    for (k, e) in arr.edges.iter().enumerate() {
        let wb = arr.below[k][0];
        let wa = wb + e.delta[0];
        if !(0..=1).contains(&wb) || !(0..=1).contains(&wa) {
            return Err(ValidityError::InvalidNesting(e.lo));
        }
    }
    // Each hole inside its own outer ring.
    for poly in polys {
        for h in &poly.holes {
            if !hole_inside(&poly.outer.0, &h.0) {
                return Err(ValidityError::InvalidNesting(h.0[0]));
            }
        }
    }
    Ok(())
}

/// `true` when the (valid, non-crossing) hole lies inside the outer ring: some hole vertex or
/// edge midpoint off the outer boundary is strictly inside.
fn hole_inside(outer: &[Point], hole: &[Point]) -> bool {
    for &p in hole {
        match locate_in_ring(outer, p) {
            Location::Inside => return true,
            Location::Outside => return false,
            Location::OnBoundary => {}
        }
    }
    // All vertices on the outer boundary: test edge midpoints in doubled coordinates.
    let dbl: Vec<Point> = outer.iter().map(|p| Point::new(2 * p.x, 2 * p.y)).collect();
    let n = hole.len();
    for i in 0..n {
        let (a, b) = (hole[i], hole[(i + 1) % n]);
        let m = Point::new(a.x + b.x, a.y + b.y);
        match locate_in_ring(&dbl, m) {
            Location::Inside => return true,
            Location::Outside => return false,
            Location::OnBoundary => {}
        }
    }
    false
}

struct UnionFind {
    parent: Vec<usize>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        UnionFind {
            parent: (0..n).collect(),
        }
    }
    fn add(&mut self) -> usize {
        self.parent.push(self.parent.len());
        self.parent.len() - 1
    }
    fn find(&mut self, mut x: usize) -> usize {
        while self.parent[x] != x {
            self.parent[x] = self.parent[self.parent[x]];
            x = self.parent[x];
        }
        x
    }
    /// Returns `false` when already connected.
    fn union(&mut self, a: usize, b: usize) -> bool {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra == rb {
            return false;
        }
        self.parent[ra] = rb;
        true
    }
}

/// Checks that a polygon set is valid **and** in the canonical form produced by this crate:
/// outer rings counter-clockwise and holes clockwise, every ring starting at its
/// lexicographically smallest vertex, holes sorted, polygons sorted by outer ring, and (when
/// `collinear_removed`) no vertex collinear with its neighbours unless it is shared with
/// another ring.
pub fn check_canonical(polys: &[Polygon], collinear_removed: bool) -> Result<(), ValidityError> {
    validate_set(polys)?;
    let mut shared: Vec<Point> = Vec::new();
    if collinear_removed {
        let mut all: Vec<Point> = polys
            .iter()
            .flat_map(|p| p.rings().flat_map(|r| r.0.iter().copied()))
            .collect();
        all.sort_unstable();
        for w in all.windows(2) {
            if w[0] == w[1] && shared.last() != Some(&w[0]) {
                shared.push(w[0]);
            }
        }
    }
    for (pi, poly) in polys.iter().enumerate() {
        for (ri, ring) in poly.rings().enumerate() {
            let id = RingId {
                polygon: pi,
                ring: ri,
            };
            if (ring_area2(&ring.0) > 0) != (ri == 0) {
                return Err(ValidityError::WrongOrientation(id));
            }
            let min = ring.0.iter().min().copied();
            if ring.0.first().copied() != min {
                return Err(ValidityError::NotCanonical(
                    "ring does not start at its minimum vertex",
                ));
            }
            if collinear_removed {
                let n = ring.0.len();
                for i in 0..n {
                    let (a, v, b) = (ring.0[(i + n - 1) % n], ring.0[i], ring.0[(i + 1) % n]);
                    if orient(a, v, b) == 0 && shared.binary_search(&v).is_err() {
                        return Err(ValidityError::NotCanonical("collinear vertex"));
                    }
                }
            }
        }
        if poly.holes.windows(2).any(|w| w[0].0 >= w[1].0) {
            return Err(ValidityError::NotCanonical("holes not sorted"));
        }
    }
    if polys.windows(2).any(|w| w[0].outer.0 >= w[1].outer.0) {
        return Err(ValidityError::NotCanonical("polygons not sorted"));
    }
    Ok(())
}
