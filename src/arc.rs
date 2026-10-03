//! Arc and circle approximation.
//!
//! Curves are converted to polylines whose distance to the true curve (the sagitta error)
//! never exceeds a tolerance, on a selectable [`Side`]:
//!
//! * [`Side::Outside`]: the approximated region contains the true region (obstacles,
//!   clearances: never under-estimates).
//! * [`Side::Inside`]: the approximated region is contained in the true region (fills).
//! * [`Side::Nearest`]: error split evenly on both sides (fewest vertices).
//!
//! The guarantee holds for the exact polyline; vertices are then rounded to the nearest
//! integer point (each coordinate moves by at most `1/2`). Vertex counts depend only on the
//! radius, the swept angle and the tolerance, and all transcendental functions come from the
//! pure-Rust `libm`, so results are bit-identical on every platform.

use crate::error::{Error, Result};
use crate::geom::{MAX_COORD, Point, Polygon, Ring, TaggedRing};
use core::f64::consts::PI;

/// Which side of the true curve the approximation may deviate to.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Side {
    /// The approximated region contains the true one.
    #[default]
    Outside,
    /// The approximated region is contained in the true one.
    Inside,
    /// Deviation on both sides, at most the tolerance each way.
    Nearest,
}

/// Arc approximation tolerance: maximum deviation (sagitta) in coordinate units and side.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ArcTol {
    /// Maximum distance between the polyline and the true arc (must be `>= 1`).
    pub tolerance: i64,
    /// Side of the deviation.
    pub side: Side,
}

impl ArcTol {
    /// Creates a tolerance.
    pub const fn new(tolerance: i64, side: Side) -> Self {
        ArcTol { tolerance, side }
    }
}

/// Maximum number of vertices generated for a single curve.
pub const MAX_ARC_VERTICES: usize = 1 << 22;

/// A circle.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Circle {
    /// Centre.
    pub center: Point,
    /// Radius (must be positive).
    pub radius: i64,
}

impl Circle {
    /// Creates a circle.
    pub const fn new(center: Point, radius: i64) -> Self {
        Circle { center, radius }
    }

    /// Approximates the circle by a canonical counter-clockwise ring (starting at its
    /// lexicographically smallest vertex). The vertex count is `n >= 4`; vertices sit at
    /// angles `2*pi*k/n` from the positive x axis (`Inside` and `Nearest`) or at the
    /// half-steps (`Outside`, so that the polygon is symmetric and tangent at the axes).
    pub fn to_ring(&self, tol: ArcTol) -> Result<Ring> {
        if self.radius <= 0 {
            return Err(Error::InvalidParameter("circle radius must be positive"));
        }
        let r = self.radius as f64;
        let step = step_angle(r, tol, true)?;
        let n = segment_count(2.0 * PI, step)?.max(4);
        let half = PI / n as f64;
        let (rv, offs) = match tol.side {
            Side::Inside => (r, 0.0),
            Side::Outside => (r / libm::cos(half), half),
            Side::Nearest => (2.0 * r / (1.0 + libm::cos(half)), 0.0),
        };
        let c = self.center;
        let mut pts = Vec::with_capacity(n);
        for k in 0..n {
            let a = offs + 2.0 * PI * k as f64 / n as f64;
            pts.push(round_pt(
                c.x as f64 + rv * libm::cos(a),
                c.y as f64 + rv * libm::sin(a),
            )?);
        }
        pts.dedup();
        if pts.len() > 1 && pts.first() == pts.last() {
            pts.pop();
        }
        let mut ring = Ring(pts);
        if let Some((i, _)) = ring.iter().enumerate().min_by_key(|(_, p)| **p) {
            ring.rotate_left(i);
        }
        Ok(ring)
    }
}

/// One element of a [`Contour`]: a straight line or an arc to a new end point.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Curve {
    /// Straight line to the point.
    Line(Point),
    /// Circular arc through `mid` to `end` (three-point form). If `end` equals the start
    /// point, this is a full circle and `mid` is the diametrically opposite point. If the
    /// three points are collinear, the arc degenerates to a line.
    Arc {
        /// A point on the arc between the start and the end.
        mid: Point,
        /// End point.
        end: Point,
    },
    /// Circular arc around `center` from the current point to `end`, counter-clockwise when
    /// `ccw`. The radius is the distance from the centre to the start point; the polyline
    /// ends exactly at `end`. If `end` equals the start point, a full circle.
    CenterArc {
        /// Centre.
        center: Point,
        /// End point.
        end: Point,
        /// Direction.
        ccw: bool,
    },
}

impl Curve {
    /// End point of the element.
    pub fn end(&self) -> Point {
        match *self {
            Curve::Line(p) => p,
            Curve::Arc { end, .. } | Curve::CenterArc { end, .. } => end,
        }
    }
}

/// A closed contour made of lines and arcs. Each element starts where the previous one
/// ends; the first element starts at the end of the last one.
pub type Contour = Vec<Curve>;

/// A region bounded by curved contours: one outer contour and holes.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Shape {
    /// Outer contour (any orientation).
    pub contour: Contour,
    /// Hole contours (any orientation).
    pub holes: Vec<Contour>,
}

impl Shape {
    /// Creates a shape.
    pub fn new(contour: Contour, holes: Vec<Contour>) -> Self {
        Shape { contour, holes }
    }

    /// Approximates the shape. The outer ring is made counter-clockwise and holes clockwise;
    /// the side of every arc is chosen relative to the material (a convex arc of the outer
    /// contour is approximated outward for `Side::Outside`, a hole's arc inward, ...).
    ///
    /// The rings are returned as computed (not passed through a boolean operation), so a
    /// self-intersecting contour stays self-intersecting; use
    /// [`union_all`](crate::union_all) to normalize.
    pub fn to_polygon(&self, tol: ArcTol) -> Result<Polygon> {
        let t = self.to_tagged(tol, &|_, _| 0)?;
        Ok(Polygon {
            outer: Ring(t[0].points.clone()),
            holes: t[1..].iter().map(|r| Ring(r.points.clone())).collect(),
        })
    }

    /// Like [`to_polygon`](Self::to_polygon) but returning tagged rings (outer first), the
    /// edges generated for element `j` of contour `i` (0 = outer, `1 + k` = hole `k`) being
    /// tagged `tag(i, j)`. Useful to reconstruct arcs after booleans and offsets.
    pub fn to_tagged(
        &self,
        tol: ArcTol,
        tag: &dyn Fn(usize, usize) -> u64,
    ) -> Result<Vec<TaggedRing>> {
        // Range-check everything first: the orientation below uses exact arithmetic.
        for c in core::iter::once(&self.contour).chain(self.holes.iter()) {
            for e in c {
                crate::error::check_point(e.end())?;
                if let Curve::Arc { mid: q, .. } | Curve::CenterArc { center: q, .. } = e {
                    crate::error::check_point(*q)?;
                }
            }
        }
        let mut out = Vec::with_capacity(1 + self.holes.len());
        for (ci, c) in core::iter::once(&self.contour)
            .chain(self.holes.iter())
            .enumerate()
        {
            let is_hole = ci > 0;
            let area = contour_area(c);
            // Material on the left when walking the contour as given?
            let ccw = area >= 0.0;
            let material_left = ccw != is_hole;
            let mut r = approx_contour(c, tol, material_left, &|j| tag(ci, j))?;
            // Canonical orientation: outer CCW, holes CW.
            let want_ccw = !is_hole;
            if ccw != want_ccw && r.points.len() > 1 {
                reverse_tagged(&mut r);
            }
            out.push(r);
        }
        Ok(out)
    }
}

/// Reverses a tagged ring's orientation, keeping each tag on its (reversed) edge.
pub(crate) fn reverse_tagged(r: &mut TaggedRing) {
    r.points.reverse();
    r.tags.reverse();
    // Edge i was points[i]->points[i+1]; after reversal edge k is points'[k]->points'[k+1]
    // which corresponds to old edge n-2-k; rotate tags by one.
    if !r.tags.is_empty() {
        r.tags.rotate_left(1);
    }
}

/// Signed area of a contour including circular segments (f64, approximate).
fn contour_area(c: &[Curve]) -> f64 {
    let Some(last) = c.last() else { return 0.0 };
    let mut cur = last.end();
    let mut a = 0.0;
    for e in c {
        let end = e.end();
        a += (cur.x as f64) * (end.y as f64) - (end.x as f64) * (cur.y as f64);
        if let Some(g) = arc_geom(cur, e) {
            // Circular segment area between chord and arc: r^2/2 (theta - sin theta),
            // signed by the sweep direction.
            let s = 0.5 * g.r * g.r * (g.sweep - libm::sin(g.sweep));
            a += 2.0 * s;
        }
        cur = end;
    }
    a / 2.0
}

#[derive(Clone, Copy, Debug)]
struct ArcGeom {
    r: f64,
    a0: f64,
    /// Signed sweep (positive = counter-clockwise).
    sweep: f64,
}

/// Geometry of an arc element starting at `start`; `None` for lines and degenerate arcs.
fn arc_geom(start: Point, e: &Curve) -> Option<ArcGeom> {
    match *e {
        Curve::Line(_) => None,
        Curve::Arc { mid, end } => {
            if end == start {
                if mid == start {
                    return None;
                }
                let cx = (start.x as f64 + mid.x as f64) / 2.0;
                let cy = (start.y as f64 + mid.y as f64) / 2.0;
                let r = libm::hypot(start.x as f64 - cx, start.y as f64 - cy);
                let a0 = libm::atan2(start.y as f64 - cy, start.x as f64 - cx);
                return Some(ArcGeom {
                    r,
                    a0,
                    sweep: 2.0 * PI,
                });
            }
            let o = crate::predicates::orient(start, mid, end);
            if o == 0 {
                return None;
            }
            // Circumcentre relative to start.
            let bx = (mid.x - start.x) as f64;
            let by = (mid.y - start.y) as f64;
            let cxr = (end.x - start.x) as f64;
            let cyr = (end.y - start.y) as f64;
            let d = 2.0 * o as f64;
            let b2 = bx * bx + by * by;
            let c2 = cxr * cxr + cyr * cyr;
            let ux = (cyr * b2 - by * c2) / d;
            let uy = (bx * c2 - cxr * b2) / d;
            let cx = start.x as f64 + ux;
            let cy = start.y as f64 + uy;
            let r = libm::hypot(ux, uy);
            let a0 = libm::atan2(-uy, -ux);
            let a1 = libm::atan2(end.y as f64 - cy, end.x as f64 - cx);
            let ccw = o > 0;
            let mut sweep = a1 - a0;
            if ccw {
                while sweep <= 0.0 {
                    sweep += 2.0 * PI;
                }
            } else {
                while sweep >= 0.0 {
                    sweep -= 2.0 * PI;
                }
            }
            Some(ArcGeom { r, a0, sweep })
        }
        Curve::CenterArc { center, end, ccw } => {
            if center == start {
                return None;
            }
            let cx = center.x as f64;
            let cy = center.y as f64;
            let r = libm::hypot(start.x as f64 - cx, start.y as f64 - cy);
            let a0 = libm::atan2(start.y as f64 - cy, start.x as f64 - cx);
            let sweep = if end == start {
                if ccw { 2.0 * PI } else { -2.0 * PI }
            } else {
                let a1 = libm::atan2(end.y as f64 - cy, end.x as f64 - cx);
                let mut s = a1 - a0;
                if ccw {
                    while s <= 0.0 {
                        s += 2.0 * PI;
                    }
                } else {
                    while s >= 0.0 {
                        s -= 2.0 * PI;
                    }
                }
                s
            };
            Some(ArcGeom { r, a0, sweep })
        }
    }
}

/// Angular step meeting the tolerance. `tangent` selects the outer (tangent-polyline)
/// construction for `Side::Outside`; for inner/outer choice by convexity see
/// [`arc_points`].
fn step_for(r: f64, tol: i64, kind: Construction) -> f64 {
    let t = tol as f64;
    // Half-step angles from numerically stable forms (no `acos` of values near 1, which
    // would round to a zero step for huge radii).
    let half = match kind {
        // cos(h/2) = 1 - t/r
        Construction::Chord if t < 2.0 * r => 2.0 * libm::asin(libm::sqrt(t / (2.0 * r))),
        // cos(h/2) = r / (r + t)
        Construction::Tangent => libm::atan(libm::sqrt(t * (2.0 * r + t)) / r),
        // cos(h/2) = (r - t) / (r + t)
        Construction::Mid if t < r => libm::atan(2.0 * libm::sqrt(r * t) / (r - t)),
        _ => PI / 4.0,
    };
    (2.0 * half).min(PI / 2.0)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Construction {
    /// Vertices on the circle (polyline inside the disk).
    Chord,
    /// Polyline tangent to the circle (outside the disk).
    Tangent,
    /// Balanced.
    Mid,
}

fn construction(side: Side, convex: bool) -> Construction {
    match (side, convex) {
        (Side::Nearest, _) => Construction::Mid,
        (Side::Outside, true) | (Side::Inside, false) => Construction::Tangent,
        (Side::Outside, false) | (Side::Inside, true) => Construction::Chord,
    }
}

fn step_angle(r: f64, tol: ArcTol, convex: bool) -> Result<f64> {
    if tol.tolerance < 1 {
        return Err(Error::InvalidParameter("arc tolerance must be >= 1"));
    }
    Ok(step_for(r, tol.tolerance, construction(tol.side, convex)))
}

fn segment_count(sweep: f64, step: f64) -> Result<usize> {
    let n = libm::ceil(sweep.abs() / step);
    if n.is_nan() || n >= MAX_ARC_VERTICES as f64 {
        return Err(Error::InvalidParameter(
            "arc tolerance too small for the radius",
        ));
    }
    Ok((n as usize).max(1))
}

#[inline]
pub(crate) fn round_pt(x: f64, y: f64) -> Result<Point> {
    let rx = libm::round(x);
    let ry = libm::round(y);
    let m = MAX_COORD as f64;
    if !(rx.abs() <= m && ry.abs() <= m) {
        let clamp = |v: f64| {
            if v.is_nan() {
                0
            } else {
                v.clamp(-(2.0 * m), 2.0 * m) as i64
            }
        };
        return Err(Error::CoordinateOutOfRange(Point::new(
            clamp(rx),
            clamp(ry),
        )));
    }
    Ok(Point::new(rx as i64, ry as i64))
}

/// Inward deviation of the chord from radius `r` (angle 0) to radius `rv` (angle `h`), for
/// arcs whose end points lie on the circle while interior vertices sit at radius `rv`.
fn end_chord_sag(r: f64, rv: f64, h: f64) -> f64 {
    let (px, py) = (r, 0.0);
    let (qx, qy) = (rv * libm::cos(h), rv * libm::sin(h));
    let (dx, dy) = (qx - px, qy - py);
    let len2 = dx * dx + dy * dy;
    // Closest point of the segment to the centre.
    let t = if len2 > 0.0 {
        (-(px * dx + py * dy) / len2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    r - libm::hypot(px + t * dx, py + t * dy)
}

/// Appends the polyline approximating an arc (excluding the start point, including `end`).
///
/// The arc starts at `(sx, sy)` (exact, possibly non-integer), has radius `r`, starts at
/// angle `a0` (as seen from the centre) and sweeps `sweep` radians. Vertices are computed
/// relative to the start point, which keeps them accurate even for huge radii. `convex`
/// tells whether the material lies on the centre's side of the arc (the region is locally
/// the disk), which decides which construction realizes `tol.side`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn arc_points(
    sx: f64,
    sy: f64,
    r: f64,
    a0: f64,
    sweep: f64,
    end: Point,
    convex: bool,
    tol: ArcTol,
    out: &mut Vec<Point>,
) -> Result<()> {
    let kind = construction(tol.side, convex);
    if tol.tolerance < 1 {
        return Err(Error::InvalidParameter("arc tolerance must be >= 1"));
    }
    let t = tol.tolerance as f64;
    let step = step_for(r, tol.tolerance, kind);
    let mut n = segment_count(sweep, step)?;
    if kind == Construction::Mid {
        // The end points stay on the circle: make sure the end chords (and a single chord)
        // also stay within the tolerance.
        loop {
            let h = (sweep / n as f64).abs();
            let rv = 2.0 * r / (1.0 + libm::cos(h / 2.0));
            let sag = if n == 1 {
                r * (1.0 - libm::cos(h / 2.0))
            } else {
                end_chord_sag(r, rv, h)
            };
            if sag <= t || n >= MAX_ARC_VERTICES {
                break;
            }
            n += 1;
        }
    }
    let h = sweep / n as f64;
    // Point at angle a0 + d on radius rv, relative to the start point:
    // (rv - r) u(a0 + d) + 2 r sin(d / 2) u_perp(a0 + d / 2).
    let at = |d: f64, rv: f64| -> Result<Point> {
        let a = a0 + d;
        let m = a0 + d / 2.0;
        let s2 = 2.0 * r * libm::sin(d / 2.0);
        let x = sx + (rv - r) * libm::cos(a) - s2 * libm::sin(m);
        let y = sy + (rv - r) * libm::sin(a) + s2 * libm::cos(m);
        round_pt(x, y)
    };
    let push = |out: &mut Vec<Point>, p: Point| {
        if out.last() != Some(&p) {
            out.push(p);
        }
    };
    match kind {
        Construction::Chord => {
            for k in 1..n {
                push(out, at(h * k as f64, r)?);
            }
        }
        Construction::Tangent => {
            let rv = r / libm::cos(h / 2.0);
            for k in 0..n {
                push(out, at(h * (k as f64 + 0.5), rv)?);
            }
        }
        Construction::Mid => {
            let rv = 2.0 * r / (1.0 + libm::cos(h / 2.0));
            for k in 1..n {
                push(out, at(h * k as f64, rv)?);
            }
        }
    }
    push(out, end);
    Ok(())
}

/// Approximates a closed contour. `material_left`: the material is on the left when walking
/// the contour in the given order.
fn approx_contour(
    c: &[Curve],
    tol: ArcTol,
    material_left: bool,
    tag: &dyn Fn(usize) -> u64,
) -> Result<TaggedRing> {
    let mut pts: Vec<Point> = Vec::new();
    let mut tags: Vec<u64> = Vec::new();
    let Some(last) = c.last() else {
        return Ok(TaggedRing::default());
    };
    let start = last.end();
    crate::error::check_point(start)?;
    let mut cur = start;
    pts.push(start);
    let mut buf = Vec::new();
    for (j, e) in c.iter().enumerate() {
        crate::error::check_point(e.end())?;
        if let Curve::Arc { mid, .. } | Curve::CenterArc { center: mid, .. } = e {
            crate::error::check_point(*mid)?;
        }
        buf.clear();
        match arc_geom(cur, e) {
            None => buf.push(e.end()),
            Some(g) => {
                // Material on the centre's side: centre is on the left of a CCW arc.
                let convex = (g.sweep > 0.0) == material_left;
                arc_points(
                    cur.x as f64,
                    cur.y as f64,
                    g.r,
                    g.a0,
                    g.sweep,
                    e.end(),
                    convex,
                    tol,
                    &mut buf,
                )?;
            }
        }
        let t = tag(j);
        for &p in &buf {
            if *pts.last().unwrap() != p {
                tags.push(t);
                pts.push(p);
            }
        }
        cur = e.end();
    }
    // Close: the last point equals the start (the last element ends there).
    if pts.len() > 1 && pts.last() == pts.first() {
        pts.pop();
        // The closing edge carries the tag of the last element.
    } else {
        // Closing edge from the last point back to the start: belongs to the last element.
        tags.push(tag(c.len() - 1));
    }
    if pts.len() == 1 {
        tags.clear();
        return Ok(TaggedRing { points: pts, tags });
    }
    // tags[i] currently labels the edge ending at pts[i + 1] (or the closing edge); this is
    // exactly edge i = pts[i] -> pts[i + 1].
    debug_assert_eq!(tags.len(), pts.len());
    Ok(TaggedRing { points: pts, tags })
}

/// Rebuilds a curved contour from a tagged ring: maximal runs of consecutive edges whose
/// tag maps to an arc (`arc_of(tag) = Some((center, ccw))`) become one
/// [`Curve::CenterArc`]; every other edge becomes a [`Curve::Line`].
///
/// Use with tags assigned per arc when approximating ([`Shape::to_tagged`]) or offsetting
/// ([`offset_shape_tagged`](crate::offset_shape_tagged)): after booleans and offsets, the
/// surviving pieces of each arc are emitted as single arcs (for example Gerber `G02`/`G03`
/// moves) instead of many short segments. The contour starts at the ring's first vertex
/// whose incoming and outgoing edges belong to different runs (or at the first vertex if
/// the whole ring is one run). The arc end points are the ring's integer vertices, so the
/// arc's radius at each end may differ from the nominal radius by the rounding (at most
/// `sqrt(2)/2`).
///
/// ```
/// use polyclip::{arcs_from_tags, ArcTol, Circle, Curve, Point, Side, TaggedRing};
/// let c = Circle::new(Point::new(0, 0), 10_000).to_ring(ArcTol::new(5, Side::Nearest)).unwrap();
/// let tagged = TaggedRing::uniform(c, 7);
/// let contour = arcs_from_tags(&tagged, &|t| (t == 7).then_some((Point::new(0, 0), true)));
/// assert_eq!(contour.len(), 1); // one full-circle arc
/// ```
pub fn arcs_from_tags(ring: &TaggedRing, arc_of: &dyn Fn(u64) -> Option<(Point, bool)>) -> Contour {
    let n = ring.points.len();
    if n == 0 || ring.tags.len() != n {
        return Vec::new();
    }
    let run_key = |i: usize| -> Option<u64> { arc_of(ring.tags[i]).map(|_| ring.tags[i]) };
    // Start at a vertex where the run changes, so no run wraps around.
    let start = (0..n).find(|&i| {
        let prev = (i + n - 1) % n;
        run_key(prev).is_none() || run_key(prev) != run_key(i)
    });
    let start = start.unwrap_or(0);
    let mut out: Contour = Vec::new();
    let mut k = 0;
    while k < n {
        let i = (start + k) % n;
        match arc_of(ring.tags[i]) {
            None => {
                out.push(Curve::Line(ring.points[(i + 1) % n]));
                k += 1;
            }
            Some((center, ccw)) => {
                let tag = ring.tags[i];
                let mut j = k;
                while j < n && ring.tags[(start + j) % n] == tag {
                    j += 1;
                }
                out.push(Curve::CenterArc {
                    center,
                    end: ring.points[(start + j) % n],
                    ccw,
                });
                k = j;
            }
        }
    }
    // As for every contour, the implicit start point is the end of the last element
    // (`points[start]`).
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::query::ring_area2;

    fn p(x: i64, y: i64) -> Point {
        Point::new(x, y)
    }

    fn max_min_radius(r: &Ring, c: Point) -> (f64, f64) {
        let mut mx: f64 = 0.0;
        let mut mn = f64::MAX;
        let n = r.len();
        for i in 0..n {
            let a = r[i];
            let b = r[(i + 1) % n];
            let da = libm::hypot((a.x - c.x) as f64, (a.y - c.y) as f64);
            mx = mx.max(da);
            // distance from centre to segment
            let (ax, ay) = ((a.x - c.x) as f64, (a.y - c.y) as f64);
            let (bx, by) = ((b.x - c.x) as f64, (b.y - c.y) as f64);
            let (dx, dy) = (bx - ax, by - ay);
            let t = (-(ax * dx + ay * dy) / (dx * dx + dy * dy)).clamp(0.0, 1.0);
            mn = mn.min(libm::hypot(ax + t * dx, ay + t * dy));
        }
        (mx, mn)
    }

    #[test]
    fn circle_sides() {
        let c = Circle::new(p(100, -50), 100_000);
        for side in [Side::Outside, Side::Inside, Side::Nearest] {
            let r = c.to_ring(ArcTol::new(10, side)).unwrap();
            assert!(ring_area2(&r) > 0);
            assert_eq!(r[0], *r.iter().min().unwrap());
            let (mx, mn) = max_min_radius(&r, c.center);
            match side {
                Side::Outside => assert!(
                    mn >= 100_000.0 - 0.75 && mx <= 100_010.0 + 0.75,
                    "{mn} {mx}"
                ),
                Side::Inside => {
                    assert!(mx <= 100_000.0 + 0.75 && mn >= 99_990.0 - 0.75, "{mn} {mx}")
                }
                Side::Nearest => assert!(mx <= 100_010.75 && mn >= 99_989.25, "{mn} {mx}"),
            }
        }
        let a = c.to_ring(ArcTol::new(10, Side::Outside)).unwrap().len();
        let b = c.to_ring(ArcTol::new(10, Side::Nearest)).unwrap().len();
        assert!(b < a);
        assert!(
            Circle::new(p(0, 0), 0)
                .to_ring(ArcTol::new(1, Side::Inside))
                .is_err()
        );
        assert!(
            Circle::new(p(0, 0), 10)
                .to_ring(ArcTol::new(0, Side::Inside))
                .is_err()
        );
        assert_eq!(
            Circle::new(p(0, 0), 1)
                .to_ring(ArcTol::new(5, Side::Inside))
                .unwrap()
                .len(),
            4
        );
    }

    #[test]
    fn shape_with_arcs() {
        // Stadium: two semicircles of radius 1000 joined by lines.
        let s = Shape::new(
            vec![
                Curve::Line(p(5000, -1000)),
                Curve::Arc {
                    mid: p(6000, 0),
                    end: p(5000, 1000),
                },
                Curve::Line(p(0, 1000)),
                Curve::Arc {
                    mid: p(-1000, 0),
                    end: p(0, -1000),
                },
            ],
            vec![vec![Curve::CenterArc {
                center: p(2500, 0),
                end: p(2700, 0),
                ccw: false,
            }]],
        );
        let tagged = s
            .to_tagged(ArcTol::new(5, Side::Outside), &|c, j| (c * 10 + j) as u64)
            .unwrap();
        assert_eq!(tagged.len(), 2);
        assert!(ring_area2(&tagged[0].points) > 0);
        assert!(ring_area2(&tagged[1].points) < 0);
        for t in &tagged {
            assert_eq!(t.points.len(), t.tags.len());
        }
        // The outer approximation contains the true stadium: area >= true area.
        let true_area = 5000.0 * 2000.0 + PI * 1000.0 * 1000.0;
        let a = ring_area2(&tagged[0].points) as f64 / 2.0;
        assert!(
            a >= true_area - 1000.0 && a <= true_area + 5.0 * 2.0 * (10000.0 + 2.0 * PI * 1000.0)
        );
        // Hole approximated inward (smaller hole) for Outside.
        let h = -(ring_area2(&tagged[1].points) as f64) / 2.0;
        assert!(h <= PI * 200.0 * 200.0 + 200.0);
        // Tags: line 0 of outer contour has tag 0, first arc tag 1, ...
        assert!(tagged[0].tags.contains(&1) && tagged[0].tags.contains(&3));
        assert!(tagged[1].tags.iter().all(|&t| t == 10));
    }
}
