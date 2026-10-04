//! Polygon and path offsetting.
//!
//! Every input ring (after normalization into canonical polygons, so outer rings are
//! counter-clockwise and holes clockwise) is turned into a raw offset curve: each edge is
//! moved by `delta` to its right (outward), and consecutive moved edges are connected at
//! each vertex either by a join (where they leave a gap) or through the original vertex
//! (where they overlap). The union of all raw curves under the [`Positive`] fill rule is the
//! offset region. This is exact in topology for any `delta`: shapes that vanish under a
//! negative offset vanish, shapes that split produce several polygons, holes that close up
//! disappear.
//!
//! Vertices of the raw curves are computed in `f64` and rounded to the nearest integer
//! point; the final union snap-rounds intersections as every boolean operation does.
//!
//! [`Positive`]: crate::FillRule::Positive

use crate::arc::Shape;
use crate::arc::{ArcTol, Circle, arc_points, round_pt};
use crate::boolean::{Boolean, FillRule, PathSource, RingSource};
use crate::error::{Error, Result};
use crate::geom::{Point, PolyTree, Polygon, PolygonSet, Ring, TaggedRing};

/// How offset edges are connected at convex corners.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Join {
    /// Circular arc around the vertex (approximated per the [`ArcTol`]).
    Round,
    /// Sharp corner where the offset edges meet, unless the corner would extend farther
    /// than `limit * |delta|` from the vertex; then it is clipped there by a line
    /// perpendicular to the corner bisector. `limit` is clamped to at least 1.
    Miter {
        /// Maximum corner distance in multiples of `|delta|`.
        limit: f64,
    },
    /// Straight line between the offset edge ends.
    Bevel,
    /// Corner clipped at distance `|delta|` from the vertex, perpendicular to the bisector.
    Square,
}

/// How the ends of open paths are shaped.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum EndCap {
    /// Half circle (tracks: stadium shapes).
    #[default]
    Round,
    /// Square extending `delta` beyond the end point.
    Square,
    /// Flat, through the end point.
    Butt,
    /// The path is closed (last vertex joined to the first) and stroked as a loop.
    Joined,
}

#[derive(Clone, Copy, Debug)]
struct V2 {
    x: f64,
    y: f64,
}

impl V2 {
    #[inline]
    fn unit(a: Point, b: Point) -> V2 {
        let dx = (b.x - a.x) as f64;
        let dy = (b.y - a.y) as f64;
        let l = libm::sqrt(dx * dx + dy * dy);
        V2 {
            x: dx / l,
            y: dy / l,
        }
    }
    /// Right-hand normal.
    #[inline]
    fn right(self) -> V2 {
        V2 {
            x: self.y,
            y: -self.x,
        }
    }
    #[inline]
    fn dot(self, o: V2) -> f64 {
        self.x * o.x + self.y * o.y
    }
    #[inline]
    fn scale(self, s: f64) -> V2 {
        V2 {
            x: self.x * s,
            y: self.y * s,
        }
    }
    #[inline]
    fn add(self, o: V2) -> V2 {
        V2 {
            x: self.x + o.x,
            y: self.y + o.y,
        }
    }
}

#[inline]
fn at(p: Point, v: V2) -> Result<Point> {
    round_pt(p.x as f64 + v.x, p.y as f64 + v.y)
}

#[inline]
fn push(out: &mut Vec<Point>, p: Point) {
    if out.last() != Some(&p) {
        out.push(p);
    }
}

/// Join used at a vertex: either the regular join or an end cap at a path end.
#[derive(Clone, Copy, Debug)]
enum VJoin {
    Join(Join),
    Cap(EndCap),
}

/// Emits the points around vertex `cur` between incoming edge `prev -> cur` and outgoing
/// edge `cur -> next`.
fn emit_vertex(
    prev: Point,
    cur: Point,
    next: Point,
    delta: f64,
    vj: VJoin,
    tol: ArcTol,
    out: &mut Vec<Point>,
) -> Result<VKind> {
    let d1 = V2::unit(prev, cur);
    let d2 = V2::unit(cur, next);
    let n1 = d1.right();
    let n2 = d2.right();
    let u1 = n1.scale(delta);
    let u2 = n2.scale(delta);
    let cr = crate::predicates::cross(
        crate::predicates::sub(cur, prev),
        crate::predicates::sub(next, cur),
    );
    let dt = crate::predicates::dot(
        crate::predicates::sub(cur, prev),
        crate::predicates::sub(next, cur),
    );
    let reversal = cr == 0 && dt < 0;
    if cr == 0 && !reversal {
        // Straight on.
        push(out, at(cur, u1)?);
        return Ok(VKind::Straight);
    }
    let convex = if reversal {
        delta > 0.0
    } else {
        (cr > 0) == (delta > 0.0)
    };
    if !convex {
        // Offset edges overlap. When they intersect well inside both (the trimmed length
        // |delta| * tan(turn / 2) is at most half of each adjacent edge, so neighbouring
        // trims cannot overlap either), use that intersection: the same region without a
        // loop. Otherwise route through the vertex itself, which is always correct.
        let c = d1.dot(d2);
        let sn = (d1.x * d2.y - d1.y * d2.x).abs();
        if c > -0.9 {
            let trim = delta.abs() * sn / (1.0 + c);
            let len1 = libm::hypot((cur.x - prev.x) as f64, (cur.y - prev.y) as f64);
            let len2 = libm::hypot((next.x - cur.x) as f64, (next.y - cur.y) as f64);
            if trim <= 0.5 * len1 && trim <= 0.5 * len2 {
                let m = n1.add(n2).scale(delta / (1.0 + c));
                push(out, at(cur, m)?);
                return Ok(VKind::Straight);
            }
        }
        push(out, at(cur, u1)?);
        push(out, cur);
        push(out, at(cur, u2)?);
        return Ok(VKind::Through);
    }
    let ad = delta.abs();
    let clip = |l: f64, out: &mut Vec<Point>| -> Result<()> {
        // Cut the corner with the line perpendicular to the bisector at distance `l`.
        let s = u1.add(u2);
        let sl = libm::sqrt(s.dot(s));
        let b = if reversal || sl < 1e-9 * ad {
            d1
        } else {
            s.scale(1.0 / sl)
        };
        let t1 = d1;
        let t2 = d2.scale(-1.0);
        let tb1 = t1.dot(b);
        let tb2 = t2.dot(b);
        if tb1 <= 1e-12 || tb2 <= 1e-12 {
            push(out, at(cur, u1)?);
            push(out, at(cur, u2)?);
            return Ok(());
        }
        let s1 = ((l - u1.dot(b)) / tb1).max(0.0);
        let s2 = ((l - u2.dot(b)) / tb2).max(0.0);
        push(out, at(cur, u1.add(t1.scale(s1)))?);
        push(out, at(cur, u2.add(t2.scale(s2)))?);
        Ok(())
    };
    let join = match vj {
        VJoin::Join(j) => j,
        VJoin::Cap(EndCap::Round) | VJoin::Cap(EndCap::Joined) => Join::Round,
        VJoin::Cap(EndCap::Square) => Join::Square,
        VJoin::Cap(EndCap::Butt) => Join::Bevel,
    };
    match join {
        Join::Bevel => {
            push(out, at(cur, u1)?);
            push(out, at(cur, u2)?);
        }
        Join::Square => clip(ad, out)?,
        Join::Miter { limit } => {
            let limit = if limit.is_nan() { 1.0 } else { limit.max(1.0) };
            let c = 1.0 + n1.dot(n2);
            // Distance of the miter point from the vertex: |delta| * sqrt(2 / (1 + cos)).
            if !reversal && c > 1e-12 && 2.0 / c <= limit * limit {
                let m = n1.add(n2).scale(delta / c);
                push(out, at(cur, m)?);
            } else {
                clip(limit * ad, out)?;
            }
        }
        Join::Round => {
            let start = at(cur, u1)?;
            let end = at(cur, u2)?;
            push(out, start);
            let a0 = libm::atan2(u1.y, u1.x);
            let mut sweep = libm::atan2(cr as f64, dt as f64);
            if reversal {
                sweep = if delta > 0.0 {
                    core::f64::consts::PI
                } else {
                    -core::f64::consts::PI
                };
            }
            let sx = cur.x as f64 + u1.x;
            let sy = cur.y as f64 + u1.y;
            // A growing offset adds a disk around the vertex (convex arc); a shrinking one
            // removes it (concave arc).
            arc_points(sx, sy, ad, a0, sweep, end, delta > 0.0, tol, out)?;
        }
    }
    Ok(VKind::Join)
}

/// How the points emitted at a vertex connect: straight on, through the vertex, or a join.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum VKind {
    Straight,
    Through,
    Join,
}

/// Appends `p` whose outgoing edge has tag `t`, merging repeated points (the later
/// outgoing tag wins).
#[inline]
fn push_t(pts: &mut Vec<Point>, tags: &mut Vec<u64>, p: Point, t: u64) {
    if pts.last() == Some(&p) {
        *tags.last_mut().unwrap() = t;
    } else {
        pts.push(p);
        tags.push(t);
    }
}

/// Raw offset curve of a closed cycle `pts` (edge `i` from `pts[i]` to `pts[i + 1]` tagged
/// `tags[i]`), with the join (or cap) at each vertex given by `vj`. Offset edges keep their
/// source tag; join edges between equally tagged edges keep that tag (an approximated arc
/// stays one arc), other joins get `corner_tag`.
#[allow(clippy::too_many_arguments)]
fn raw_cycle(
    pts: &[Point],
    tags: &[u64],
    delta: f64,
    vj: impl Fn(usize) -> VJoin,
    tol: ArcTol,
    corner_tag: u64,
    out: &mut Vec<TaggedRing>,
) -> Result<()> {
    let n = pts.len();
    if n < 2 {
        return Ok(());
    }
    let mut rp: Vec<Point> = Vec::with_capacity(n * 2);
    let mut rt: Vec<u64> = Vec::with_capacity(n * 2);
    let mut vbuf: Vec<Point> = Vec::new();
    for i in 0..n {
        let prev = pts[(i + n - 1) % n];
        let next = pts[(i + 1) % n];
        let (tin, tout) = (tags[(i + n - 1) % n], tags[i]);
        vbuf.clear();
        let kind = emit_vertex(prev, pts[i], next, delta, vj(i), tol, &mut vbuf)?;
        let k = vbuf.len();
        for (j, &p) in vbuf.iter().enumerate() {
            let t = if j + 1 == k {
                tout
            } else {
                match kind {
                    VKind::Through => {
                        if j == 0 {
                            tin
                        } else {
                            tout
                        }
                    }
                    _ if tin == tout => tin,
                    _ => corner_tag,
                }
            };
            push_t(&mut rp, &mut rt, p, t);
        }
    }
    while rp.len() > 1 && rp.first() == rp.last() {
        rp.pop();
        rt.pop();
    }
    if rp.len() >= 3 {
        out.push(TaggedRing {
            points: rp,
            tags: rt,
        });
    }
    Ok(())
}

/// Removes consecutive duplicates (and a closing duplicate when `closed`), keeping for each
/// kept point the tag of its outgoing edge.
fn dedup_tagged(pts: &[Point], tags: &[u64], closed: bool) -> (Vec<Point>, Vec<u64>) {
    let mut v: Vec<Point> = Vec::with_capacity(pts.len());
    let mut t: Vec<u64> = Vec::with_capacity(pts.len());
    for (i, &p) in pts.iter().enumerate() {
        let tag = tags.get(i).copied().unwrap_or(0);
        if v.last() == Some(&p) {
            *t.last_mut().unwrap() = tag;
        } else {
            v.push(p);
            t.push(tag);
        }
    }
    if closed {
        while v.len() > 1 && v.first() == v.last() {
            v.pop();
            t.pop();
        }
    }
    (v, t)
}

fn check_delta(delta: i64) -> Result<()> {
    if delta.unsigned_abs() > 2 * crate::MAX_COORD as u64 {
        return Err(Error::InvalidParameter("offset delta too large"));
    }
    Ok(())
}

/// Offsets a region by `delta`, propagating edge tags, and returns the nesting tree.
///
/// Like [`offset_tree`], but every output edge carries the tag of the input edge it was
/// offset from. Join edges between two edges with the same tag (for example consecutive
/// segments of an approximated arc, see [`Shape::to_tagged`](crate::Shape::to_tagged)) keep
/// that tag, so an offset arc is still recognizable as one arc; joins at corners between
/// differently tagged edges are tagged `corner_tag`. Intersections created by the final
/// union keep the tags of the edges they lie on.
pub fn offset_tagged(
    input: &(impl RingSource + ?Sized),
    delta: i64,
    join: Join,
    tol: ArcTol,
    corner_tag: u64,
) -> Result<PolyTree> {
    check_delta(delta)?;
    // Already canonical polygons (e.g. the output of an earlier boolean) are a fixed point
    // of normalization: offset them directly.
    let mut polys: Vec<&Polygon> = Vec::new();
    if delta != 0 && input.visit_polygons(&mut |p| polys.push(p)) && looks_canonical(&polys) {
        let mut rings: Vec<TaggedRing> = Vec::new();
        let mut zeros: Vec<u64> = Vec::new();
        for p in polys {
            for r in p.rings() {
                zeros.clear();
                zeros.resize(r.len(), 0);
                raw_cycle(
                    &r.0,
                    &zeros,
                    delta as f64,
                    |_| VJoin::Join(join),
                    tol,
                    corner_tag,
                    &mut rings,
                )?;
            }
        }
        return Boolean::new()
            .subject(&rings, FillRule::Positive)
            .execute_tree();
    }
    let norm = Boolean::new()
        .subject(input, FillRule::NonZero)
        .execute_tree()?;
    if delta == 0 {
        return Ok(norm);
    }
    offset_canonical_tree(&norm, delta, join, tol, corner_tag)
}

/// Offsets a canonical tree (the output of a boolean) without re-normalizing it.
fn offset_canonical_tree(
    norm: &PolyTree,
    delta: i64,
    join: Join,
    tol: ArcTol,
    corner_tag: u64,
) -> Result<PolyTree> {
    let mut rings: Vec<TaggedRing> = Vec::new();
    for n in &norm.nodes {
        raw_cycle(
            &n.ring,
            &n.tags,
            delta as f64,
            |_| VJoin::Join(join),
            tol,
            corner_tag,
            &mut rings,
        )?;
    }
    Boolean::new()
        .subject(&rings, FillRule::Positive)
        .execute_tree()
}

/// Offsets a region by `delta` (positive grows, negative shrinks), returning the full
/// nesting tree.
///
/// The input rings are first normalized with the non-zero fill rule (so any orientation and
/// overlap is accepted). Holes are handled naturally: they shrink when the region grows and
/// grow when it shrinks. A negative offset larger than a feature's half-width removes it;
/// necks narrower than `2 * |delta|` split the shape. `delta == 0` just normalizes.
///
/// Polygons already in canonical form (such as the output of a boolean or another offset)
/// are recognized in linear time and not normalized again; such input is assumed to be
/// valid (see [`validate`](crate::validate)).
pub fn offset_tree(
    input: &(impl RingSource + ?Sized),
    delta: i64,
    join: Join,
    tol: ArcTol,
) -> Result<PolyTree> {
    offset_tagged(input, delta, join, tol, 0)
}

/// Offsets a region by `delta` (positive grows, negative shrinks). See [`offset_tree`].
///
/// ```
/// use polyclip::{offset, ArcTol, Join, Ring, Side};
/// let sq = Ring::from([(0, 0), (100, 0), (100, 100), (0, 100)]);
/// let grown = offset(&sq, 10, Join::Miter { limit: 2.0 }, ArcTol::new(1, Side::Outside)).unwrap();
/// assert_eq!(grown[0].outer, Ring::from([(-10, -10), (110, -10), (110, 110), (-10, 110)]));
/// let gone = offset(&sq, -50, Join::Round, ArcTol::new(1, Side::Inside)).unwrap();
/// assert!(gone.is_empty());
/// ```
pub fn offset(
    input: &(impl RingSource + ?Sized),
    delta: i64,
    join: Join,
    tol: ArcTol,
) -> Result<PolygonSet> {
    Ok(offset_tree(input, delta, join, tol)?.to_polygon_set())
}

/// Offsets a curved [`Shape`] by `delta`.
///
/// The shape is approximated with half the tolerance and the result offset with the other
/// half, both on `tol.side`, so the total deviation from the true offset shape stays within
/// `tol.tolerance` on the requested side (up to integer rounding). Typical use: pad shapes
/// inflated by a clearance with `Side::Outside`, so the obstacle never under-estimates.
///
/// ```
/// use polyclip::{offset_shape, ArcTol, Circle, Curve, Join, Point, Shape, Side};
/// let pad = Shape::new(vec![Curve::CenterArc { center: Point::new(0, 0), end: Point::new(500, 0), ccw: true }]
///     .into_iter().chain([Curve::Line(Point::new(500, 0))]).collect(), vec![]);
/// let obstacle = offset_shape(&pad, 200, Join::Round, ArcTol::new(10, Side::Outside)).unwrap();
/// assert_eq!(obstacle.len(), 1);
/// ```
pub fn offset_shape(shape: &Shape, delta: i64, join: Join, tol: ArcTol) -> Result<PolygonSet> {
    Ok(offset_shape_tagged(shape, delta, join, tol, &|_, _| 0, 0)?.to_polygon_set())
}

/// Like [`offset_shape`], with edge tags: edges coming from element `j` of contour `i` (0 =
/// outer, `1 + k` = hole `k`) are tagged `tag(i, j)`, corner joins `corner_tag`.
pub fn offset_shape_tagged(
    shape: &Shape,
    delta: i64,
    join: Join,
    tol: ArcTol,
    tag: &dyn Fn(usize, usize) -> u64,
    corner_tag: u64,
) -> Result<PolyTree> {
    let half = ArcTol::new((tol.tolerance / 2).max(1), tol.side);
    let rings = shape.to_tagged(half, tag)?;
    offset_tagged(&rings, delta, join, half, corner_tag)
}

/// Morphological opening: shrink by `d`, then grow by `d`, with round joins. Removes every
/// part narrower than `2 * d` (minimum-width enforcement for copper zones) while leaving
/// wide parts essentially unchanged (convex corners get rounded with radius `d`).
pub fn opening(input: &(impl RingSource + ?Sized), d: i64, tol: ArcTol) -> Result<PolygonSet> {
    let d = d
        .checked_abs()
        .ok_or(Error::InvalidParameter("offset delta too large"))?;
    // The intermediate tree is canonical: grow it directly.
    let shrunk = offset_tagged(input, -d, Join::Round, tol, 0)?;
    Ok(offset_canonical_tree(&shrunk, d, Join::Round, tol, 0)?.to_polygon_set())
}

/// Morphological closing: grow by `d`, then shrink by `d`, with round joins. Fills gaps and
/// notches narrower than `2 * d`.
pub fn closing(input: &(impl RingSource + ?Sized), d: i64, tol: ArcTol) -> Result<PolygonSet> {
    let d = d
        .checked_abs()
        .ok_or(Error::InvalidParameter("offset delta too large"))?;
    let grown = offset_tagged(input, d, Join::Round, tol, 0)?;
    Ok(offset_canonical_tree(&grown, -d, Join::Round, tol, 0)?.to_polygon_set())
}

/// Cheap (linear) check that polygons are in the canonical form produced by booleans:
/// rings with at least three vertices starting at their smallest one, no repeated or
/// collinear vertices, outer rings counter-clockwise, holes clockwise, holes and polygons
/// strictly sorted, coordinates in range. Normalizing such input returns it unchanged
/// (validity is assumed, not checked: the offset of crossing input is still valid output).
fn looks_canonical(polys: &[&Polygon]) -> bool {
    let ring_ok = |r: &Ring, outer: bool| -> bool {
        let p = &r.0;
        let n = p.len();
        if n < 3 || !p.iter().all(|q| q.in_range()) || p[1..].iter().any(|q| *q <= p[0]) {
            return false;
        }
        if (0..n).any(|i| crate::predicates::orient(p[i], p[(i + 1) % n], p[(i + 2) % n]) == 0) {
            return false;
        }
        (crate::query::ring_area2(p) > 0) == outer
    };
    polys.windows(2).all(|w| w[0].outer.0 < w[1].outer.0)
        && polys.iter().all(|p| {
            ring_ok(&p.outer, true)
                && p.holes.iter().all(|h| ring_ok(h, false))
                && p.holes.windows(2).all(|w| w[0].0 < w[1].0)
        })
}

/// Offsets open paths with edge tags (see [`offset_tagged`] for the tagging rules; end caps
/// are tagged `corner_tag`). Returns the nesting tree.
pub fn offset_paths_tagged(
    paths: &(impl PathSource + ?Sized),
    delta: i64,
    join: Join,
    cap: EndCap,
    tol: ArcTol,
    corner_tag: u64,
) -> Result<PolyTree> {
    if delta < 0 {
        return Err(Error::InvalidParameter("path offset must be non-negative"));
    }
    check_delta(delta)?;
    let mut rings: Vec<TaggedRing> = Vec::new();
    let mut err: Option<Error> = None;
    paths.visit_paths(&mut |pts, tags| {
        if err.is_some() {
            return;
        }
        if let Some(&p) = pts.iter().find(|p| !p.in_range()) {
            err = Some(Error::CoordinateOutOfRange(p));
            return;
        }
        if delta == 0 {
            return;
        }
        let tags = tags.unwrap_or(&[]);
        if let Err(e) = raw_path(
            pts,
            tags,
            delta as f64,
            join,
            cap,
            tol,
            corner_tag,
            &mut rings,
        ) {
            err = Some(e);
        }
    });
    if let Some(e) = err {
        return Err(e);
    }
    Boolean::new()
        .subject(&rings, FillRule::Positive)
        .execute_tree()
}

/// Offsets open paths by `delta >= 0` on both sides (strokes them with width `2 * delta`),
/// with `join` at interior vertices and `cap` at both ends. Returns the nesting tree.
///
/// A path with a single distinct point becomes a circle (round cap), a square (square cap)
/// or nothing (butt cap). With [`EndCap::Joined`] the path is treated as closed.
pub fn offset_paths_tree(
    paths: &(impl PathSource + ?Sized),
    delta: i64,
    join: Join,
    cap: EndCap,
    tol: ArcTol,
) -> Result<PolyTree> {
    offset_paths_tagged(paths, delta, join, cap, tol, 0)
}

/// Offsets open paths. See [`offset_paths_tree`].
///
/// ```
/// use polyclip::{offset_paths, ArcTol, EndCap, Join, Path, Side};
/// let track = Path::from([(0, 0), (1000, 0)]);
/// let s = offset_paths(&track, 100, Join::Round, EndCap::Butt, ArcTol::new(1, Side::Outside)).unwrap();
/// assert_eq!(s[0].outer.signed_area2(), 2 * 1000 * 200);
/// ```
pub fn offset_paths(
    paths: &(impl PathSource + ?Sized),
    delta: i64,
    join: Join,
    cap: EndCap,
    tol: ArcTol,
) -> Result<PolygonSet> {
    Ok(offset_paths_tree(paths, delta, join, cap, tol)?.to_polygon_set())
}

#[allow(clippy::too_many_arguments)]
fn raw_path(
    pts: &[Point],
    tags: &[u64],
    delta: f64,
    join: Join,
    cap: EndCap,
    tol: ArcTol,
    corner_tag: u64,
    rings: &mut Vec<TaggedRing>,
) -> Result<()> {
    let closed = cap == EndCap::Joined;
    // For a closed path the closing edge carries the last tag.
    let (v, t) = dedup_tagged(pts, tags, closed);
    match v.len() {
        0 => return Ok(()),
        1 => {
            let p = v[0];
            let d = delta as i64;
            let ring = match cap {
                EndCap::Round | EndCap::Joined => Circle::new(p, d).to_ring(tol)?,
                EndCap::Square => {
                    for q in [Point::new(p.x - d, p.y - d), Point::new(p.x + d, p.y + d)] {
                        crate::error::check_point(q)?;
                    }
                    Ring::from([
                        (p.x - d, p.y - d),
                        (p.x + d, p.y - d),
                        (p.x + d, p.y + d),
                        (p.x - d, p.y + d),
                    ])
                }
                EndCap::Butt => return Ok(()),
            };
            rings.push(TaggedRing::uniform(ring, corner_tag));
            return Ok(());
        }
        _ => {}
    }
    if closed && v.len() >= 3 {
        // Both sides of the loop: the ring and its reverse, each offset outward.
        raw_cycle(&v, &t, delta, |_| VJoin::Join(join), tol, corner_tag, rings)?;
        let mut rev = TaggedRing { points: v, tags: t };
        crate::arc::reverse_tagged(&mut rev);
        raw_cycle(
            &rev.points,
            &rev.tags,
            delta,
            |_| VJoin::Join(join),
            tol,
            corner_tag,
            rings,
        )?;
        return Ok(());
    }
    // Walk forward then back: p0 .. pn .. p1, with caps at the two turnarounds.
    let n = v.len();
    let mut ring: Vec<Point> = Vec::with_capacity(2 * n);
    ring.extend_from_slice(&v);
    ring.extend(v[1..n - 1].iter().rev());
    let m = ring.len();
    // Edge i < n - 1 is path edge i; edge i >= n - 1 is path edge 2n - 3 - i reversed.
    let rtags: Vec<u64> = (0..m)
        .map(|i| if i < n - 1 { t[i] } else { t[2 * n - 3 - i] })
        .collect();
    let cap = if closed { EndCap::Round } else { cap };
    raw_cycle(
        &ring,
        &rtags,
        delta,
        |i| {
            if i == 0 || i == n - 1 {
                VJoin::Cap(cap)
            } else {
                VJoin::Join(join)
            }
        },
        tol,
        corner_tag,
        rings,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arc::Side;
    use crate::geom::Path;
    use crate::validate::check_canonical;

    fn sq(x0: i64, y0: i64, x1: i64, y1: i64) -> Ring {
        Ring::from([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])
    }
    fn area(ps: &PolygonSet) -> f64 {
        ps.iter().map(|p| p.signed_area2() as f64 / 2.0).sum()
    }
    const TOL: ArcTol = ArcTol::new(2, Side::Nearest);

    #[test]
    fn square_joins() {
        let s = sq(0, 0, 1000, 1000);
        let m = offset(&s, 100, Join::Miter { limit: 2.0 }, TOL).unwrap();
        assert_eq!(m, vec![crate::Polygon::from(sq(-100, -100, 1100, 1100))]);
        let b = offset(&s, 100, Join::Bevel, TOL).unwrap();
        assert_eq!(area(&b), 1200.0 * 1200.0 - 4.0 * 100.0 * 100.0 / 2.0);
        let q = offset(&s, 100, Join::Square, TOL).unwrap();
        // Clipped at distance 100 along the diagonal: cut corners.
        let a = area(&q);
        assert!(a > area(&b) && a < 1200.0 * 1200.0, "{a}");
        let r = offset(&s, 100, Join::Round, TOL).unwrap();
        let ar = area(&r);
        let exact = 1000.0 * 1000.0 + 4.0 * 1000.0 * 100.0 + core::f64::consts::PI * 100.0 * 100.0;
        assert!((ar - exact).abs() < 0.01 * exact, "{ar} {exact}");
        for x in [&m, &b, &q, &r] {
            assert_eq!(check_canonical(x, true), Ok(()));
        }
        // Miter limit 1 behaves like square.
        assert_eq!(offset(&s, 100, Join::Miter { limit: 1.0 }, TOL).unwrap(), q);
    }

    #[test]
    fn shrink_split_vanish() {
        // Dumbbell: two 100x100 squares joined by a 10-wide neck.
        let d = crate::union_all(
            &vec![
                sq(0, 0, 100, 100),
                sq(200, 0, 300, 100),
                sq(100, 45, 200, 55),
            ],
            FillRule::NonZero,
        )
        .unwrap();
        let s = offset(&d, -10, Join::Miter { limit: 2.0 }, TOL).unwrap();
        assert_eq!(
            s,
            vec![
                crate::Polygon::from(sq(10, 10, 90, 90)),
                crate::Polygon::from(sq(210, 10, 290, 90))
            ]
        );
        assert!(offset(&d, -60, Join::Round, TOL).unwrap().is_empty());
        let o = opening(&d, 6, TOL).unwrap();
        assert_eq!(o.len(), 2);
    }

    #[test]
    fn holes() {
        let p = crate::Polygon::new(
            sq(0, 0, 100, 100),
            vec![Ring::from([(40, 40), (40, 60), (60, 60), (60, 40)])],
        );
        let g = offset(&p, 5, Join::Miter { limit: 2.0 }, TOL).unwrap();
        assert_eq!(g[0].outer, sq(-5, -5, 105, 105));
        assert_eq!(
            g[0].holes,
            vec![Ring::from([(45, 45), (45, 55), (55, 55), (55, 45)])]
        );
        let g = offset(&p, 10, Join::Miter { limit: 2.0 }, TOL).unwrap();
        assert!(g[0].holes.is_empty());
        let s = offset(&p, -5, Join::Miter { limit: 2.0 }, TOL).unwrap();
        assert_eq!(s[0].outer, sq(5, 5, 95, 95));
        assert_eq!(
            s[0].holes,
            vec![Ring::from([(35, 35), (35, 65), (65, 65), (65, 35)])]
        );
    }

    #[test]
    fn paths() {
        let t = Path::from([(0, 0), (1000, 0)]);
        let round = offset_paths(&t, 100, Join::Round, EndCap::Round, TOL).unwrap();
        let exact = 1000.0 * 200.0 + core::f64::consts::PI * 100.0 * 100.0;
        assert!((area(&round) - exact).abs() < 0.01 * exact);
        let sq_cap = offset_paths(&t, 100, Join::Round, EndCap::Square, TOL).unwrap();
        assert_eq!(sq_cap[0].outer, sq(-100, -100, 1100, 100));
        let l = Path::from([(0, 0), (1000, 0), (1000, 1000)]);
        let lm = offset_paths(&l, 100, Join::Miter { limit: 2.0 }, EndCap::Butt, TOL).unwrap();
        assert_eq!(area(&lm), 1100.0 * 200.0 + 900.0 * 200.0);
        let dot = offset_paths(
            &Path::from([(5, 5), (5, 5)]),
            100,
            Join::Round,
            EndCap::Round,
            TOL,
        )
        .unwrap();
        assert_eq!(dot.len(), 1);
        assert!(
            offset_paths(&Path::from([(5, 5)]), 100, Join::Round, EndCap::Butt, TOL)
                .unwrap()
                .is_empty()
        );
        let lp = Path::from([(0, 0), (100, 0), (100, 100), (0, 100)]);
        let j = offset_paths(&lp, 10, Join::Miter { limit: 2.0 }, EndCap::Joined, TOL).unwrap();
        assert_eq!(j[0].outer, sq(-10, -10, 110, 110));
        assert_eq!(
            j[0].holes,
            vec![Ring::from([(10, 10), (10, 90), (90, 90), (90, 10)])]
        );
        // Self-overlapping zig-zag still gives a valid result.
        let z = Path::from([(0, 0), (100, 0), (0, 5), (100, 10)]);
        let zz = offset_paths(&z, 20, Join::Round, EndCap::Round, TOL).unwrap();
        assert_eq!(check_canonical(&zz, true), Ok(()));
    }

    #[test]
    fn concave_round() {
        // L shape grown with round joins: reflex corner stays sharp.
        let l = Ring::from([
            (0, 0),
            (200, 0),
            (200, 100),
            (100, 100),
            (100, 200),
            (0, 200),
        ]);
        let g = offset(&l, 10, Join::Round, TOL).unwrap();
        assert!(g[0].outer.contains(&Point::new(110, 110)));
        let s = offset(&l, -10, Join::Round, TOL).unwrap();
        assert!(s[0].outer.contains(&Point::new(10, 10)));
        assert_eq!(check_canonical(&s, true), Ok(()));
    }
}

#[cfg(test)]
mod tag_tests {
    use super::*;
    use crate::arc::{Curve, Side};

    #[test]
    fn tags_through_offset() {
        let p = Point::new;
        // Stadium: line, arc, line, arc.
        let s = Shape::new(
            vec![
                Curve::Line(p(50_000, -10_000)),
                Curve::Arc {
                    mid: p(60_000, 0),
                    end: p(50_000, 10_000),
                },
                Curve::Line(p(0, 10_000)),
                Curve::Arc {
                    mid: p(-10_000, 0),
                    end: p(0, -10_000),
                },
            ],
            vec![],
        );
        let t = offset_shape_tagged(
            &s,
            2_000,
            Join::Round,
            ArcTol::new(20, Side::Outside),
            &|_, j| 100 + j as u64,
            999,
        )
        .unwrap();
        assert_eq!(t.nodes.len(), 1);
        let n = &t.nodes[0];
        assert_eq!(n.tags.len(), n.ring.len());
        // Every source element survives as a tag; arcs dominate the vertex count.
        for tag in [100, 101, 102, 103] {
            assert!(n.tags.contains(&tag), "missing {tag}: {:?}", n.tags);
        }
        // Edges tagged with an arc lie (approximately) on the offset circle of radius 12000.
        let centre = |tag| {
            if tag == 101 {
                (50_000.0, 0.0)
            } else {
                (0.0, 0.0)
            }
        };
        for (i, &tag) in n.tags.iter().enumerate() {
            if tag == 101 || tag == 103 {
                let q = n.ring[i];
                let (cx, cy) = centre(tag);
                let r = ((q.x as f64 - cx).powi(2) + (q.y as f64 - cy).powi(2)).sqrt();
                assert!((r - 12_000.0).abs() < 30.0, "tag {tag} vertex {q:?} r {r}");
            }
        }
        // Same region as the untagged version (which also merges collinear edges whose tags
        // differ, so the vertex lists may differ).
        let untagged =
            offset_shape(&s, 2_000, Join::Round, ArcTol::new(20, Side::Outside)).unwrap();
        assert_eq!(untagged[0].outer.signed_area2(), n.ring.signed_area2());
    }
}
