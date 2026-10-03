//! Arc-preserving boolean operations on curved [`Shape`]s.
//!
//! See [`curved_boolean`] for the algorithm and its guarantees.

use crate::arc::{ArcTol, Contour, Curve, Shape, Side, round_pt};
use crate::boolean::{Boolean, FillRule, Op};
use crate::error::{Error, Result};
use crate::geom::{Point, Polygon, TaggedPolygon, TaggedRing};
use crate::predicates::{dist2, orient};
use crate::query::ring_area2;
use crate::validate::{ValidityError, validate_set};
use core::f64::consts::{PI, TAU};
use std::collections::{HashMap, HashSet};

/// Safety margin (in coordinate units) added to every band-separation test.
const MARGIN: f64 = 2.0;
/// Half-width of the band around a straight input edge (snap rounding moves edges by at
/// most `sqrt(2)/2`).
const LINE_W: f64 = 2.0;
/// Minimum sine of the crossing angle for a crossing to be reconstructed exactly.
const MIN_SIN: f64 = 1.0 / 32.0;
/// Radius above which an arc is not reconstructed (kept as its approximation).
const MAX_ARC_RADIUS: f64 = (1u64 << 42) as f64;
/// Maximum distance between a centre-form arc's end point and its circle for the arc to be
/// reconstructed.
const MAX_END_OFF: f64 = 2.0;
/// Number of verification rounds before falling back to the plain polygon result.
const MAX_ROUNDS: usize = 16;

/// Boolean operation on curved shapes that keeps arcs as arcs.
///
/// Computes `op(subject, clip)` where each operand is the union (under `rule`) of its
/// shapes, and returns the result as [`Shape`]s whose curved parts are real arcs
/// ([`Curve::CenterArc`]) lying on the circles of the input arcs, so that outputs such as
/// Gerber `G02`/`G03` or IPC-2581 arcs can be written directly. Straight parts are
/// [`Curve::Line`]s.
///
/// # What "native" means here
///
/// The intersection point of two circles (or of a circle and a line) is irrational, so an
/// exact curved result with integer vertices does not exist in general. This function
/// computes the *topology* with the crate's exact polygon boolean on sided arc
/// approximations, and then rebuilds the curves:
///
/// 1. Every input element gets its own tag; arcs are approximated with tolerance
///    `tol.tolerance` on the side that makes the polygon result err towards `tol.side`
///    (for [`Op::Difference`] with [`Side::Inside`], for example, the subject is
///    approximated inward and the clip outward). The exact polygon boolean of these
///    approximations gives the reference result `P`.
/// 2. Around every input curve lies a *band* of half-width `tol.tolerance + 2` (`2` for
///    lines) containing both the true curve and every polygon edge derived from it. Where
///    the band of an arc is separated from all other bands (by `2` more units), the true
///    boundary there is that arc, and the polygon run of `P` is replaced by the true arc.
/// 3. Where two curves cross transversally (crossing angle with sine at least `1/32`) with
///    exactly one true crossing in the region where their bands overlap, and no third band
///    reaches that region, the vertex of `P` is moved to the true crossing point (rounded
///    to the nearest integer point) and both curves run up to it. A junction between two
///    consecutive elements of an input contour is kept at the input vertex. (If an arc
///    cannot leave the moved vertex, because its band is not clean right after it, the
///    move is undone: a straight edge from the true crossing could cut through the arc.)
/// 4. Everywhere else (near-tangencies, several crossings close together, crossings close
///    to a third curve, coincident curves of different circles) the output follows `P`
///    exactly; the arc resumes, with a short radial connector, as soon as its band is
///    clean again.
/// 5. The result is verified: each rebuilt ring's curved area must match `P` within the
///    band bound, and the output approximated with `ArcTol::new(tol.tolerance,
///    Side::Nearest)` must pass [`validate_set`] (full circles rebuilt in a clean band are
///    valid by construction and skipped). Rings and polygons failing a check are replaced
///    by their (always valid) `P` counterpart, made of lines only.
///
/// # Guarantees
///
/// Let `t = tol.tolerance`, `R` the exact result of the operation on the true curves and
/// `C` the returned shapes.
///
/// * **Arcs**: every output arc is a [`Curve::CenterArc`] whose centre is the centre of an
///   input arc (exactly when that centre is an integer point, which it always is for
///   [`Curve::CenterArc`] input; three-point arcs with a non-integral circumcentre get the
///   nearest integer point). Its start and end points lie within one unit of that input
///   circle (`sqrt(2)/2` for integral centres), so its radius is the input radius up to
///   rounding. Arcs never sweep more than half a turn, except full circles, which are a
///   single element ending at its own start point.
/// * **Deviation**: every point where `C` and `R` differ lies within `t + 2` of the
///   boundary of an input shape. (A Hausdorff bound alone cannot be given: like any
///   approximation, a sliver of `R` thinner than `t` may be absent from `C`.)
/// * **Side**: with [`Side::Inside`], `C` is contained in `R` up to rounding: points of `C`
///   outside `R` are within 2 units of the boundary of an input shape. With
///   [`Side::Outside`], `C` contains `R` in the same sense. For example
///   `zone − obstacles` with `Side::Inside` never intrudes into a true obstacle by more than
///   rounding, and approximations that touch or nearly touch are resolved conservatively.
///   [`Side::Nearest`] gives no side guarantee. The side guarantee assumes the shapes of an
///   operand combine by union ([`FillRule::NonZero`] or [`FillRule::Positive`]); with
///   other rules it holds for each shape but not where shapes of one operand overlap.
/// * **Validity**: `P` is valid, and every output polygon that keeps arcs passed the
///   checks of step 5, so `C` approximated with `ArcTol::new(t, Side::Nearest)` is always a
///   valid polygon set (the property tests also check tolerance 1). Rebuilt arcs keep at
///   least `4` units from every other curve of `C` except where they meet it at a vertex,
///   and moved crossings have a crossing angle with sine `>= 1/32`. As for any curved
///   shape, a coarse *one-sided* approximation can self-intersect at a sharp corner (two
///   curves meeting tangentially, for instance); [`union_all`](crate::union_all)
///   normalizes it.
/// * **Which arcs survive**: an input arc survives wherever its band is clean, up to the
///   rounded true crossings of step 3; a circle with a clean band all around (an isolated
///   via or round pad) comes out as one full-circle element starting at the input's own
///   start point.
/// * **Degenerate cases**: arcs on the same exact circle (same integer centre and squared
///   radius) are merged, so coincident arcs come out as one arc; tangent arcs, arcs
///   touching lines and crossings at grazing angles keep the conservative polygon geometry
///   locally. Arc pieces whose sagitta is below half a unit are written as lines. Arcs with
///   a radius above `2^42`, and centre-form arcs whose end point is more than 2 units off
///   their circle, are not reconstructed (their approximation is kept).
/// * **Determinism**: the result depends only on the input (no hashing order, no
///   threads), and is bit-identical across platforms. Shapes come in the canonical order
///   of the polygon result; outer contours are counter-clockwise and holes clockwise; each
///   contour starts at its lexicographically smallest vertex.
///
/// Shapes of each operand are approximated with [`Shape::to_tagged`], so outer contours
/// count counter-clockwise and holes clockwise whatever their input orientation. The fill
/// rule applies to the oriented rings of each operand. [`Op::Xor`] with `Side::Inside` or
/// `Side::Outside` is computed as `(A − B) ∪ (B − A)`, each difference with its own sides.
///
/// # Errors
///
/// [`Error::InvalidParameter`] when `tol.tolerance < 1` or the tolerance is too small for a
/// radius, [`Error::CoordinateOutOfRange`] for out-of-range input (or approximations
/// leaving the range), [`Error::TooLarge`] for huge inputs.
///
/// # Performance
///
/// The cost is that of the polygon boolean on the approximations, plus the band tests
/// (a sweep over bounding boxes), the rebuild and the validation of the rebuilt output
/// (isolated full circles excluded). A 100 mm board with rounded corners minus 5 000 round
/// pads (1 um tolerance) takes about 100 ms, against about 58 ms for the polygon boolean
/// alone (`cargo bench -- curved_zone`).
///
/// ```
/// use polyclip::*;
/// let p = Point::new;
/// // A 10 x 10 mm board with a 1 mm corner radius, minus a round via.
/// let board = Shape::new(vec![
///     Curve::Line(p(9_000_000, 0)),
///     Curve::CenterArc { center: p(9_000_000, 1_000_000), end: p(10_000_000, 1_000_000), ccw: true },
///     Curve::Line(p(10_000_000, 10_000_000)),
///     Curve::Line(p(0, 10_000_000)),
///     Curve::Line(p(0, 0)),
/// ], vec![]);
/// let via = Shape::new(vec![Curve::CenterArc { center: p(5_000_000, 5_000_000), end: p(5_300_000, 5_000_000), ccw: true }], vec![]);
/// let fill = curved_boolean(Op::Difference, &[board], &[via], FillRule::NonZero,
///     ArcTol::new(1_000, Side::Inside)).unwrap();
/// assert_eq!(fill.len(), 1);
/// // The corner arc and the via survive as arcs on their own circles.
/// let arcs = |c: &Contour| c.iter().filter(|e| matches!(e, Curve::CenterArc { .. })).count();
/// assert_eq!(arcs(&fill[0].contour), 1);
/// assert_eq!(fill[0].holes.len(), 1);
/// assert_eq!(fill[0].holes[0], vec![Curve::CenterArc { center: p(5_000_000, 5_000_000), end: p(5_300_000, 5_000_000), ccw: false }]);
/// ```
pub fn curved_boolean(
    op: Op,
    subject: &[Shape],
    clip: &[Shape],
    rule: FillRule,
    tol: ArcTol,
) -> Result<Vec<Shape>> {
    if tol.tolerance < 1 {
        return Err(Error::InvalidParameter("arc tolerance must be >= 1"));
    }
    if tol.tolerance > crate::geom::MAX_COORD {
        return Err(Error::InvalidParameter("arc tolerance too large"));
    }
    let mut ctx = Ctx::new(tol.tolerance as f64);
    let sub_bases = ctx.add_operand(subject);
    let clip_bases = ctx.add_operand(clip);

    // Polygon result on sided approximations.
    let side = tol.side;
    let approx = |shapes: &[Shape], bases: &[Vec<usize>], s: Side| -> Result<Vec<TaggedRing>> {
        let at = ArcTol::new(tol.tolerance, s);
        let mut out = Vec::new();
        for (shape, b) in shapes.iter().zip(bases) {
            out.extend(shape.to_tagged(at, &|i, j| (b[i] + j + 1) as u64)?);
        }
        Ok(out)
    };
    let mut all_rings: Vec<Vec<TaggedRing>> = Vec::new();
    let p: Vec<TaggedPolygon> = if op == Op::Xor && side != Side::Nearest {
        let a_s = approx(subject, &sub_bases, side)?;
        let a_f = approx(subject, &sub_bases, flip(side))?;
        let b_s = approx(clip, &clip_bases, side)?;
        let b_f = approx(clip, &clip_bases, flip(side))?;
        let d1 = Boolean::new()
            .subject(&a_s, rule)
            .clip(&b_f, rule)
            .op(Op::Difference)
            .execute_tagged()?;
        let d2 = Boolean::new()
            .subject(&b_s, rule)
            .clip(&a_f, rule)
            .op(Op::Difference)
            .execute_tagged()?;
        all_rings.extend([a_s, a_f, b_s, b_f]);
        Boolean::new()
            .subject(&d1, FillRule::NonZero)
            .clip(&d2, FillRule::NonZero)
            .op(Op::Union)
            .execute_tagged()?
    } else {
        let (ss, cs) = match op {
            Op::Difference => (side, flip(side)),
            _ => (side, side),
        };
        let a = approx(subject, &sub_bases, ss)?;
        let b = approx(clip, &clip_bases, cs)?;
        let r = Boolean::new()
            .subject(&a, rule)
            .clip(&b, rule)
            .op(op)
            .execute_tagged()?;
        all_rings.extend([a, b]);
        r
    };
    if p.is_empty() {
        return Ok(Vec::new());
    }
    ctx.add_opaque_blockers(&all_rings);
    drop(all_rings);
    ctx.compute_blocked();
    Ok(ctx.reconstruct(&p))
}

fn flip(s: Side) -> Side {
    match s {
        Side::Inside => Side::Outside,
        Side::Outside => Side::Inside,
        Side::Nearest => Side::Nearest,
    }
}

// ---------------------------------------------------------------------------------------
// Small f64 vector helpers.

type V = [f64; 2];

#[inline]
fn fv(p: Point) -> V {
    [p.x as f64, p.y as f64]
}
#[inline]
fn vsub(a: V, b: V) -> V {
    [a[0] - b[0], a[1] - b[1]]
}
#[inline]
fn vadd(a: V, b: V) -> V {
    [a[0] + b[0], a[1] + b[1]]
}
#[inline]
fn vmul(a: V, s: f64) -> V {
    [a[0] * s, a[1] * s]
}
#[inline]
fn vdot(a: V, b: V) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
#[inline]
fn vcross(a: V, b: V) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}
#[inline]
fn vlen(a: V) -> f64 {
    libm::hypot(a[0], a[1])
}
#[inline]
fn ang(a: V) -> f64 {
    libm::atan2(a[1], a[0])
}
#[inline]
fn unit(a: f64) -> V {
    [libm::cos(a), libm::sin(a)]
}
/// Normalizes an angle to `[0, 2 pi)`.
fn norm_tau(a: f64) -> f64 {
    let r = a - TAU * libm::floor(a / TAU);
    if (0.0..TAU).contains(&r) { r } else { 0.0 }
}
/// Normalizes an angle to `(-pi, pi]`.
fn norm_pi(a: f64) -> f64 {
    let r = norm_tau(a);
    if r > PI { r - TAU } else { r }
}

// ---------------------------------------------------------------------------------------
// Input elements.

/// Geometry of one input element, mirroring how [`Shape::to_tagged`] approximates it.
#[derive(Clone, Copy, Debug)]
enum Raw {
    /// Straight edge (possibly zero-length).
    Line(Point, Point),
    /// Circular arc that can be reconstructed.
    Arc {
        c: V,
        r: f64,
        a0: f64,
        /// Signed sweep (positive = counter-clockwise).
        sweep: f64,
        start: Point,
        end: Point,
        /// Exact circle (integer centre, squared radius) when the centre is integral.
        key: Option<(Point, i128)>,
        c_out: Point,
    },
    /// Arc kept as its approximation.
    Opaque,
}

fn sweep_to(a0: f64, a1: f64, ccw: bool) -> f64 {
    let mut s = a1 - a0;
    if ccw {
        while s <= 0.0 {
            s += TAU;
        }
    } else {
        while s >= 0.0 {
            s -= TAU;
        }
    }
    s
}

fn classify(start: Point, e: &Curve) -> Raw {
    match *e {
        Curve::Line(p) => Raw::Line(start, p),
        Curve::Arc { mid, end } => {
            if end == start {
                if mid == start {
                    return Raw::Line(start, end);
                }
                let sx2 = start.x as i128 + mid.x as i128;
                let sy2 = start.y as i128 + mid.y as i128;
                let c = [sx2 as f64 / 2.0, sy2 as f64 / 2.0];
                let r = vlen(vsub(fv(start), c));
                let a0 = ang(vsub(fv(start), c));
                let key = (sx2 % 2 == 0 && sy2 % 2 == 0).then(|| {
                    let cp = Point::new((sx2 / 2) as i64, (sy2 / 2) as i64);
                    (cp, dist2(start, cp))
                });
                return arc_raw(c, r, a0, TAU, start, end, key);
            }
            let o = orient(start, mid, end);
            if o == 0 {
                return Raw::Line(start, end);
            }
            let bx = (mid.x - start.x) as i128;
            let by = (mid.y - start.y) as i128;
            let cx = (end.x - start.x) as i128;
            let cy = (end.y - start.y) as i128;
            let d = 2 * o;
            let b2 = bx * bx + by * by;
            let c2 = cx * cx + cy * cy;
            let uxn = cy * b2 - by * c2;
            let uyn = bx * c2 - cx * b2;
            let ux = uxn as f64 / d as f64;
            let uy = uyn as f64 / d as f64;
            let c = [start.x as f64 + ux, start.y as f64 + uy];
            let r = libm::hypot(ux, uy);
            let a0 = libm::atan2(-uy, -ux);
            let a1 = libm::atan2(end.y as f64 - c[1], end.x as f64 - c[0]);
            let sweep = sweep_to(a0, a1, o > 0);
            let key = if uxn % d == 0 && uyn % d == 0 {
                let (qx, qy) = (uxn / d, uyn / d);
                let lim = 1i128 << 42;
                (qx.abs() <= lim && qy.abs() <= lim).then(|| {
                    let cp = Point::new(start.x + qx as i64, start.y + qy as i64);
                    (cp, dist2(start, cp))
                })
            } else {
                None
            };
            arc_raw(c, r, a0, sweep, start, end, key)
        }
        Curve::CenterArc { center, end, ccw } => {
            if center == start {
                return Raw::Line(start, end);
            }
            let c = fv(center);
            let r = vlen(vsub(fv(start), c));
            let a0 = ang(vsub(fv(start), c));
            let sweep = if end == start {
                if ccw { TAU } else { -TAU }
            } else {
                sweep_to(a0, ang(vsub(fv(end), c)), ccw)
            };
            if (vlen(vsub(fv(end), c)) - r).abs() > MAX_END_OFF {
                return Raw::Opaque;
            }
            arc_raw(
                c,
                r,
                a0,
                sweep,
                start,
                end,
                Some((center, dist2(start, center))),
            )
        }
    }
}

fn arc_raw(
    c: V,
    r: f64,
    a0: f64,
    sweep: f64,
    start: Point,
    end: Point,
    key: Option<(Point, i128)>,
) -> Raw {
    if !r.is_finite() || r <= 0.0 || r > MAX_ARC_RADIUS || !c[0].is_finite() || !c[1].is_finite() {
        return Raw::Opaque;
    }
    let c_out = match key {
        Some((p, _)) => p,
        None => match round_pt(c[0], c[1]) {
            Ok(p) => p,
            Err(_) => return Raw::Opaque,
        },
    };
    if !c_out.in_range() {
        return Raw::Opaque;
    }
    Raw::Arc {
        c,
        r,
        a0,
        sweep,
        start,
        end,
        key,
        c_out,
    }
}

/// A circle that arcs are reconstructed on.
#[derive(Clone, Copy, Debug)]
struct Circ {
    c: V,
    r: f64,
    c_out: Point,
    /// An integer point exactly on the circle (start of a full-circle input element).
    full_start: Option<Point>,
}

/// Geometry used for band separation tests.
#[derive(Clone, Copy, Debug)]
enum Prim {
    Seg {
        a: V,
        b: V,
    },
    Arc {
        c: V,
        r: f64,
        a0: f64,
        sweep: f64,
        e0: V,
        e1: V,
    },
}

#[derive(Clone, Copy, Debug)]
struct Blocker {
    curve: u32,
    prim: Prim,
    w: f64,
}

/// Parametrized curve on which blocked parameter intervals are computed.
#[derive(Clone, Copy, Debug)]
enum Tgt {
    /// Full circle, parameter = angle in `[0, 2 pi)`.
    Circle { c: V, r: f64 },
    /// Segment from `a` along unit direction `d`, parameter = distance in `[0, len]`.
    Seg { a: V, d: V, len: f64 },
}

impl Tgt {
    fn periodic(&self) -> bool {
        matches!(self, Tgt::Circle { .. })
    }
    fn at(&self, s: f64) -> V {
        match *self {
            Tgt::Circle { c, r } => vadd(c, vmul(unit(s), r)),
            Tgt::Seg { a, d, .. } => vadd(a, vmul(d, s)),
        }
    }
    fn param(&self, p: V) -> f64 {
        match *self {
            Tgt::Circle { c, .. } => norm_tau(ang(vsub(p, c))),
            Tgt::Seg { a, d, .. } => vdot(vsub(p, a), d),
        }
    }
    /// Parameter tolerance corresponding to about half a unit of length.
    fn eps(&self) -> f64 {
        match *self {
            Tgt::Circle { r, .. } => (0.5 / r.max(0.5)).min(0.01),
            Tgt::Seg { .. } => 0.5,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Kind {
    Circle(usize),
    Line,
    Opaque,
}

#[derive(Clone, Debug)]
struct CurveInfo {
    kind: Kind,
    tgt: Option<Tgt>,
    w: f64,
    /// Parameter intervals where another curve's band comes too close, with that curve.
    blocked: Vec<(f64, f64, u32)>,
}

/// A junction moved to the true crossing of its two curves.
#[derive(Clone, Copy, Debug)]
struct Snap {
    y: Point,
    t: V,
    comp_u: (f64, f64),
    comp_v: (f64, f64),
}

/// End of a run: where it starts or stops, and which blocked component is resolved there.
#[derive(Clone, Copy, Debug)]
struct End {
    p: Point,
    /// Moved away from the vertex of `P` (to the true crossing).
    moved: bool,
    t: Option<V>,
    exempt: Option<(u32, (f64, f64))>,
}

/// Why a run could not be rebuilt: a moved junction at its start or end, or anything else.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Fail {
    Start,
    End,
    Ring,
}

fn blame(s: &End, e: &End) -> Fail {
    if s.moved {
        Fail::Start
    } else if e.moved {
        Fail::End
    } else {
        Fail::Ring
    }
}

enum RingOut {
    Curved(Contour),
    Drop,
}

struct Ctx {
    tol: f64,
    w_arc: f64,
    /// Curve of each element (index = tag - 1).
    elem_curve: Vec<u32>,
    elem_raw: Vec<Raw>,
    curves: Vec<CurveInfo>,
    circles: Vec<Circ>,
    circle_keys: HashMap<(Point, i128), usize>,
    /// Curve id of each circle.
    circle_curve: Vec<u32>,
    blockers: Vec<Blocker>,
    /// `(vertex, min curve, max curve)` for every vertex between consecutive input elements.
    originals: HashSet<(Point, u32, u32)>,
}

impl Ctx {
    fn new(tol: f64) -> Self {
        Ctx {
            tol,
            w_arc: tol + 2.0,
            elem_curve: Vec::new(),
            elem_raw: Vec::new(),
            curves: Vec::new(),
            circles: Vec::new(),
            circle_keys: HashMap::new(),
            circle_curve: Vec::new(),
            blockers: Vec::new(),
            originals: HashSet::new(),
        }
    }

    fn new_curve(&mut self, kind: Kind, tgt: Option<Tgt>, w: f64) -> u32 {
        self.curves.push(CurveInfo {
            kind,
            tgt,
            w,
            blocked: Vec::new(),
        });
        (self.curves.len() - 1) as u32
    }

    /// Registers the elements of an operand; returns per shape the element base index of
    /// each contour.
    fn add_operand(&mut self, shapes: &[Shape]) -> Vec<Vec<usize>> {
        let mut bases = Vec::with_capacity(shapes.len());
        for s in shapes {
            let mut b = Vec::with_capacity(1 + s.holes.len());
            for c in core::iter::once(&s.contour).chain(s.holes.iter()) {
                b.push(self.elem_curve.len());
                self.add_contour(c);
            }
            bases.push(b);
        }
        bases
    }

    fn add_contour(&mut self, c: &[Curve]) {
        let Some(last) = c.last() else { return };
        let mut cur = last.end();
        let first = self.elem_curve.len();
        for e in c {
            // Out-of-range input is rejected later by the approximation; avoid overflow in
            // the exact arithmetic here.
            let ok = cur.in_range()
                && e.end().in_range()
                && match e {
                    Curve::Arc { mid: q, .. } | Curve::CenterArc { center: q, .. } => q.in_range(),
                    Curve::Line(_) => true,
                };
            let raw = if ok {
                classify(cur, e)
            } else {
                Raw::Line(cur, cur)
            };
            let cid = match raw {
                Raw::Line(a, b) => {
                    let tgt = (a != b).then(|| {
                        let (fa, fb) = (fv(a), fv(b));
                        let len = vlen(vsub(fb, fa));
                        Tgt::Seg {
                            a: fa,
                            d: vmul(vsub(fb, fa), 1.0 / len),
                            len,
                        }
                    });
                    let id = self.new_curve(Kind::Line, tgt, LINE_W);
                    if a != b {
                        self.blockers.push(Blocker {
                            curve: id,
                            prim: Prim::Seg { a: fv(a), b: fv(b) },
                            w: LINE_W,
                        });
                    }
                    id
                }
                Raw::Opaque => self.new_curve(Kind::Opaque, None, self.w_arc),
                Raw::Arc {
                    c,
                    r,
                    a0,
                    sweep,
                    start,
                    end,
                    key,
                    c_out,
                } => {
                    let full = sweep.abs() >= TAU;
                    let circ = match key.and_then(|k| self.circle_keys.get(&k).copied()) {
                        Some(i) => i,
                        None => {
                            let (cc, rr) = match key {
                                Some((p, r2)) => (fv(p), libm::sqrt(r2 as f64)),
                                None => (c, r),
                            };
                            self.circles.push(Circ {
                                c: cc,
                                r: rr,
                                c_out,
                                full_start: None,
                            });
                            let i = self.circles.len() - 1;
                            if let Some(k) = key {
                                self.circle_keys.insert(k, i);
                            }
                            let id = self.new_curve(
                                Kind::Circle(i),
                                Some(Tgt::Circle { c: cc, r: rr }),
                                self.w_arc,
                            );
                            self.circle_curve.push(id);
                            i
                        }
                    };
                    if full && key.is_some() && self.circles[circ].full_start.is_none() {
                        self.circles[circ].full_start = Some(start);
                    }
                    let id = self.circle_curve[circ];
                    self.blockers.push(Blocker {
                        curve: id,
                        prim: Prim::Arc {
                            c,
                            r,
                            a0,
                            sweep,
                            e0: fv(start),
                            e1: fv(end),
                        },
                        w: self.w_arc,
                    });
                    id
                }
            };
            self.elem_curve.push(cid);
            self.elem_raw.push(raw);
            cur = e.end();
        }
        let n = c.len();
        if n >= 2 {
            for (j, e) in c.iter().enumerate() {
                let (a, b) = (first + j, first + (j + 1) % n);
                let (ca, cb) = (self.elem_curve[a], self.elem_curve[b]);
                if ca != cb {
                    let v = e.end();
                    self.originals.insert((v, ca.min(cb), ca.max(cb)));
                }
            }
        }
    }

    /// Opaque elements block with the segments of their approximation.
    fn add_opaque_blockers(&mut self, rings: &[Vec<TaggedRing>]) {
        for set in rings {
            for r in set {
                let n = r.points.len();
                for i in 0..n {
                    let Some(&tag) = r.tags.get(i) else { continue };
                    let e = (tag as usize).wrapping_sub(1);
                    if e < self.elem_raw.len() && matches!(self.elem_raw[e], Raw::Opaque) {
                        let (a, b) = (r.points[i], r.points[(i + 1) % n]);
                        if a != b {
                            self.blockers.push(Blocker {
                                curve: self.elem_curve[e],
                                prim: Prim::Seg { a: fv(a), b: fv(b) },
                                w: self.w_arc,
                            });
                        }
                    }
                }
            }
        }
    }

    fn curve_of_tag(&self, tag: u64) -> Option<u32> {
        let e = (tag as usize).checked_sub(1)?;
        self.elem_curve.get(e).copied()
    }

    /// Computes, for every curve, the parameter intervals where another curve's band comes
    /// within the separation distance.
    fn compute_blocked(&mut self) {
        #[derive(Clone, Copy)]
        struct Item {
            bb: [f64; 4],
            idx: usize,
            target: bool,
        }
        let mut items: Vec<Item> = Vec::new();
        for (i, cv) in self.curves.iter().enumerate() {
            let Some(t) = cv.tgt else { continue };
            let e = cv.w + MARGIN;
            let bb = match t {
                Tgt::Circle { c, r } => [c[0] - r - e, c[1] - r - e, c[0] + r + e, c[1] + r + e],
                Tgt::Seg { a, d, len } => {
                    let b = vadd(a, vmul(d, len));
                    [
                        a[0].min(b[0]) - e,
                        a[1].min(b[1]) - e,
                        a[0].max(b[0]) + e,
                        a[1].max(b[1]) + e,
                    ]
                }
            };
            items.push(Item {
                bb,
                idx: i,
                target: true,
            });
        }
        for (i, b) in self.blockers.iter().enumerate() {
            let e = b.w;
            let bb = prim_bbox(&b.prim);
            items.push(Item {
                bb: [bb[0] - e, bb[1] - e, bb[2] + e, bb[3] + e],
                idx: i,
                target: false,
            });
        }
        items.sort_by(|a, b| {
            a.bb[0]
                .total_cmp(&b.bb[0])
                .then(a.target.cmp(&b.target))
                .then(a.idx.cmp(&b.idx))
        });
        let mut act_t: Vec<Item> = Vec::new();
        let mut act_b: Vec<Item> = Vec::new();
        let mut pairs: Vec<(usize, usize)> = Vec::new();
        for it in &items {
            let x = it.bb[0];
            act_t.retain(|a| a.bb[2] >= x);
            act_b.retain(|a| a.bb[2] >= x);
            let yo = |a: &Item| a.bb[1] <= it.bb[3] && it.bb[1] <= a.bb[3];
            if it.target {
                for b in act_b.iter().filter(|b| yo(b)) {
                    pairs.push((it.idx, b.idx));
                }
                act_t.push(*it);
            } else {
                for t in act_t.iter().filter(|t| yo(t)) {
                    pairs.push((t.idx, it.idx));
                }
                act_b.push(*it);
            }
        }
        pairs.sort_unstable();
        for (ti, bi) in pairs {
            let b = self.blockers[bi];
            if b.curve as usize == ti {
                continue;
            }
            let cv = &self.curves[ti];
            let Some(t) = cv.tgt else { continue };
            let wsum = cv.w + b.w + MARGIN;
            for (lo, hi) in blocked_intervals(&t, &b.prim, wsum) {
                self.curves[ti].blocked.push((lo, hi, b.curve));
            }
        }
        for cv in &mut self.curves {
            cv.blocked.sort_by(|a, b| {
                a.0.total_cmp(&b.0)
                    .then(a.1.total_cmp(&b.1))
                    .then(a.2.cmp(&b.2))
            });
        }
    }

    /// Rebuilds the curved result from the polygon result `p`.
    fn reconstruct(&self, p: &[TaggedPolygon]) -> Vec<Shape> {
        // Vertices shared between rings (or repeated within one) are never moved, unless
        // every edge at the vertex comes from one circle (coincident approximations of the
        // same circle, whose true boundary is that circle alone).
        const MIXED: u32 = u32::MAX;
        let mut count: HashMap<Point, (u32, u32)> = HashMap::new();
        for poly in p {
            for r in core::iter::once(&poly.outer).chain(poly.holes.iter()) {
                let n = r.points.len();
                for (i, &q) in r.points.iter().enumerate() {
                    let mut kind = MIXED;
                    if r.tags.len() == n
                        && let (Some(a), Some(b)) = (
                            self.curve_of_tag(r.tags[(i + n - 1) % n]),
                            self.curve_of_tag(r.tags[i]),
                        )
                        && a == b
                        && matches!(self.curves[a as usize].kind, Kind::Circle(_))
                    {
                        kind = a;
                    }
                    let e = count.entry(q).or_insert((0, kind));
                    e.0 += 1;
                    if e.1 != kind {
                        e.1 = MIXED;
                    }
                }
            }
        }
        let hard = |q: Point| count.get(&q).is_some_and(|&(k, c)| k > 1 && c == MIXED);

        // Junction decisions, with uniqueness per crossing region.
        let mut rings: Vec<&TaggedRing> = Vec::new();
        for poly in p {
            rings.push(&poly.outer);
            rings.extend(poly.holes.iter());
        }
        type Key = (u32, u32, u64, u64);
        let mut cands: Vec<((usize, usize), Snap, Key)> = Vec::new();
        let mut groups: HashMap<Key, (Point, bool)> = HashMap::new();
        for (ri, r) in rings.iter().enumerate() {
            let n = r.points.len();
            if n < 3 || r.tags.len() != n {
                continue;
            }
            for i in 0..n {
                let (Some(u), Some(v)) = (
                    self.curve_of_tag(r.tags[(i + n - 1) % n]),
                    self.curve_of_tag(r.tags[i]),
                ) else {
                    continue;
                };
                if u == v {
                    continue;
                }
                let x = r.points[i];
                if let Some((s, key)) = self.decide(x, u, v) {
                    // Several vertices of `P` for one crossing, or one shared with another
                    // ring: keep them all.
                    let g = groups.entry(key).or_insert((x, false));
                    if g.0 != x || hard(x) {
                        g.1 = true;
                    }
                    cands.push(((ri, i), s, key));
                }
            }
        }
        let mut snaps: HashMap<(usize, usize), Snap> = HashMap::new();
        for (at, s, key) in cands {
            if !groups[&key].1 {
                snaps.insert(at, s);
            }
        }

        // Build each ring.
        let mut out: Vec<Option<Shape>> = Vec::with_capacity(p.len());
        let mut ri = 0;
        for poly in p {
            let mut rs: Vec<Option<RingOut>> = Vec::new();
            for r in core::iter::once(&poly.outer).chain(poly.holes.iter()) {
                rs.push(self.build_ring(ri, r, &snaps, &hard));
                ri += 1;
            }
            let ring_contour = |k: usize, o: Option<RingOut>| -> Option<Contour> {
                let tr = if k == 0 {
                    &poly.outer
                } else {
                    &poly.holes[k - 1]
                };
                match o {
                    Some(RingOut::Curved(c)) => {
                        if self.area_ok(&c, tr) {
                            Some(c)
                        } else {
                            Some(lines(&tr.points))
                        }
                    }
                    Some(RingOut::Drop) => None,
                    None => Some(lines(&tr.points)),
                }
            };
            let mut it = rs.into_iter().enumerate();
            let (_, o) = it.next().unwrap_or((0, None));
            match ring_contour(0, o) {
                None => out.push(None),
                Some(outer) => {
                    let holes: Vec<Contour> = it.filter_map(|(k, o)| ring_contour(k, o)).collect();
                    // An outer ring and a hole rebuilt to the same full circle bound an
                    // empty region (two approximations of one circle).
                    let same_circle = |h: &Contour| match (outer.as_slice(), h.as_slice()) {
                        (
                            [
                                Curve::CenterArc {
                                    center: c1,
                                    end: e1,
                                    ..
                                },
                            ],
                            [
                                Curve::CenterArc {
                                    center: c2,
                                    end: e2,
                                    ..
                                },
                            ],
                        ) => c1 == c2 && e1 == e2,
                        _ => false,
                    };
                    if holes.iter().any(same_circle) {
                        out.push(None);
                    } else {
                        out.push(Some(Shape::new(outer, holes)));
                    }
                }
            }
        }

        // Verification: approximate and validate; replace offending polygons by `P`.
        let fallback = |i: usize| -> Shape {
            Shape::new(
                lines(&p[i].outer.points),
                p[i].holes.iter().map(|h| lines(&h.points)).collect(),
            )
        };
        let vt = ArcTol::new(self.tol as i64, Side::Nearest);
        let mut reverted = vec![false; p.len()];
        let mut validated = false;
        for round in 0..=MAX_ROUNDS {
            let mut polys: Vec<Polygon> = Vec::new();
            let mut owner: Vec<usize> = Vec::new();
            let mut retry = false;
            for (i, s) in out.iter().enumerate() {
                let Some(s) = s else { continue };
                if reverted[i] {
                    polys.push(p[i].clone().into_polygon());
                } else {
                    // Full circles rebuilt in a clean band cannot meet anything else: leave
                    // them out of the (costly) check.
                    let holes: Vec<Contour> = s
                        .holes
                        .iter()
                        .filter(|h| !self.certain(h))
                        .cloned()
                        .collect();
                    if holes.is_empty() && self.certain(&s.contour) {
                        continue;
                    }
                    match Shape::new(s.contour.clone(), holes).to_polygon(vt) {
                        Ok(pg) => polys.push(pg),
                        Err(_) => {
                            reverted[i] = true;
                            retry = true;
                            polys.push(p[i].clone().into_polygon());
                        }
                    }
                }
                owner.push(i);
            }
            if retry {
                continue;
            }
            match validate_set(&polys) {
                Ok(()) => {
                    validated = true;
                    break;
                }
                Err(e) => {
                    let hits = error_polygons(&e, &polys);
                    let mut changed = false;
                    for k in hits {
                        let i = owner[k];
                        if !reverted[i] {
                            reverted[i] = true;
                            changed = true;
                        }
                    }
                    if !changed || round == MAX_ROUNDS {
                        break;
                    }
                }
            }
        }
        if !validated {
            // Give up on arcs entirely: `P` is valid.
            reverted.iter_mut().for_each(|r| *r = true);
            out = (0..p.len()).map(|i| Some(fallback(i))).collect();
        }
        out.into_iter()
            .enumerate()
            .filter_map(|(i, s)| {
                let s = if reverted[i] { Some(fallback(i)) } else { s }?;
                Some(Shape::new(
                    canonical_start(s.contour),
                    s.holes.into_iter().map(canonical_start).collect(),
                ))
            })
            .collect()
    }

    /// Decides whether the junction `x` between curves `u` (incoming) and `v` (outgoing) can
    /// be moved to the true crossing.
    fn decide(&self, x: Point, u: u32, v: u32) -> Option<(Snap, (u32, u32, u64, u64))> {
        let (cu, cv) = (&self.curves[u as usize], &self.curves[v as usize]);
        if cu.kind == Kind::Opaque || cv.kind == Kind::Opaque {
            return None;
        }
        if cu.kind == Kind::Line && cv.kind == Kind::Line {
            return None;
        }
        let (tu, tv) = (cu.tgt?, cv.tgt?);
        let xf = fv(x);
        let comp_u = component(&cu.blocked, v, tu.param(xf), &tu)?;
        let comp_v = component(&cv.blocked, u, tv.param(xf), &tv)?;
        if others_overlap(&cu.blocked, v, comp_u, &tu)
            || others_overlap(&cv.blocked, u, comp_v, &tv)
        {
            return None;
        }
        let original = self.originals.contains(&(x, u.min(v), u.max(v)));
        let (t, y) = if original {
            (xf, x)
        } else {
            let mut found = None;
            let mut count = 0;
            for (q, sin) in intersections(&tu, &tv) {
                let pu = shift_into(tu.param(q), comp_u, &tu);
                let pv = shift_into(tv.param(q), comp_v, &tv);
                if inside(pu, comp_u, tu.eps()) && inside(pv, comp_v, tv.eps()) {
                    count += 1;
                    found = Some((q, sin));
                }
            }
            let (q, sin) = found?;
            if count != 1 || sin < MIN_SIN {
                return None;
            }
            (q, round_pt(q[0], q[1]).ok()?)
        };
        let key_comp = if u < v { comp_u } else { comp_v };
        Some((
            Snap {
                y,
                t,
                comp_u,
                comp_v,
            },
            (
                u.min(v),
                u.max(v),
                key_comp.0.to_bits(),
                key_comp.1.to_bits(),
            ),
        ))
    }

    fn build_ring(
        &self,
        ri: usize,
        r: &TaggedRing,
        snaps: &HashMap<(usize, usize), Snap>,
        hard: &dyn Fn(Point) -> bool,
    ) -> Option<RingOut> {
        let n = r.points.len();
        if n < 3 || r.tags.len() != n {
            return None;
        }
        let cids: Vec<u32> = r
            .tags
            .iter()
            .map(|&t| self.curve_of_tag(t))
            .collect::<Option<Vec<_>>>()?;
        let pts = &r.points;
        let junctions: Vec<usize> = (0..n)
            .filter(|&i| cids[(i + n - 1) % n] != cids[i])
            .collect();
        let mut out: Contour = Vec::new();
        if junctions.is_empty() {
            let cid = cids[0];
            let Kind::Circle(ci) = self.curves[cid as usize].kind else {
                return None;
            };
            let circ = &self.circles[ci];
            let cv = &self.curves[cid as usize];
            let any_hard = pts.iter().any(|&q| hard(q));
            // Total turn around the centre.
            let mut total = 0.0;
            let mut lo = f64::MAX;
            let mut hi = f64::MIN;
            let mut a = ang(vsub(fv(pts[0]), circ.c));
            let base = a;
            for k in 1..=n {
                let b = ang(vsub(fv(pts[k % n]), circ.c));
                total += norm_pi(b - a);
                lo = lo.min(base + total);
                hi = hi.max(base + total);
                a = b;
            }
            if total.abs() > PI {
                if !any_hard && cv.blocked.is_empty() {
                    let s = match circ.full_start {
                        Some(s) => s,
                        None => round_pt(circ.c[0] + circ.r, circ.c[1]).ok()?,
                    };
                    out.push(Curve::CenterArc {
                        center: circ.c_out,
                        end: s,
                        ccw: total > 0.0,
                    });
                    return Some(RingOut::Curved(out));
                }
                // Partly blocked: start at a vertex where the circle is not clean.
                let tgt = cv.tgt?;
                let start = (0..n).find(|&i| {
                    hard(pts[i]) || {
                        let th = tgt.param(fv(pts[i]));
                        !clean(&cv.blocked, th, th, &tgt, None, None)
                    }
                });
                let start = start.or_else(|| {
                    (0..n).find(|&i| {
                        let a = tgt.param(fv(pts[i]));
                        let b = a + norm_pi(tgt.param(fv(pts[(i + 1) % n])) - a);
                        !clean(&cv.blocked, a.min(b), a.max(b), &tgt, None, None)
                    })
                })?;
                let mut run: Vec<Point> = (0..=n).map(|k| pts[(start + k) % n]).collect();
                run.dedup();
                let e = End {
                    p: pts[start],
                    moved: false,
                    t: None,
                    exempt: None,
                };
                self.emit_run(&run, cid, e, e, hard, &mut out).ok()?;
                return Some(RingOut::Curved(out));
            }
            // A sliver between approximations of one circle: nothing in the true result.
            let tgt = cv.tgt?;
            if !any_hard && clean(&cv.blocked, lo, hi, &tgt, None, None) {
                return Some(RingOut::Drop);
            }
            return None;
        }
        let m = junctions.len();
        // Junctions whose move turned out unsafe for this ring are kept at `P`'s vertex.
        let mut disabled = vec![false; m];
        let mut done = false;
        let end_at = |j: usize, incoming: bool, disabled: &[bool]| -> End {
            let i = junctions[j % m];
            match snaps.get(&(ri, i)).filter(|_| !disabled[j % m]) {
                Some(s) => {
                    let (other, comp) = if incoming {
                        (cids[i], s.comp_u)
                    } else {
                        (cids[(i + n - 1) % n], s.comp_v)
                    };
                    End {
                        p: s.y,
                        moved: s.y != pts[i],
                        t: Some(s.t),
                        exempt: Some((other, comp)),
                    }
                }
                None => End {
                    p: pts[i],
                    moved: false,
                    t: None,
                    exempt: None,
                },
            }
        };
        'attempt: for _ in 0..=m {
            out.clear();
            for j in 0..m {
                let (a, b) = (junctions[j], junctions[(j + 1) % m]);
                let len = if b > a { b - a } else { b + n - a };
                let run: Vec<Point> = (0..=len).map(|k| pts[(a + k) % n]).collect();
                let s = end_at(j, false, &disabled);
                let e = end_at(j + 1, true, &disabled);
                match self.emit_run(&run, cids[a], s, e, hard, &mut out) {
                    Ok(()) => {}
                    Err(Fail::Start) => {
                        disabled[j] = true;
                        continue 'attempt;
                    }
                    Err(Fail::End) => {
                        disabled[(j + 1) % m] = true;
                        continue 'attempt;
                    }
                    Err(Fail::Ring) => return None,
                }
            }
            done = true;
            break;
        }
        if !done || out.is_empty() {
            return None;
        }
        Some(RingOut::Curved(out))
    }

    /// Appends the elements of one run (all edges from one curve), from `s` (already the
    /// current point) to `e`.
    fn emit_run(
        &self,
        run: &[Point],
        cid: u32,
        s: End,
        e: End,
        hard: &dyn Fn(Point) -> bool,
        out: &mut Contour,
    ) -> core::result::Result<(), Fail> {
        let mut cur = s.p;
        let push_line = |out: &mut Contour, cur: &mut Point, q: Point| {
            if *cur != q {
                out.push(Curve::Line(q));
                *cur = q;
            }
        };
        let cv = &self.curves[cid as usize];
        if run.len() < 2 {
            return Err(Fail::Ring);
        }
        let m = run.len() - 1;
        match cv.kind {
            Kind::Opaque => {
                if s.t.is_some() || e.t.is_some() {
                    return Err(blame(&s, &e));
                }
                for &q in &run[1..] {
                    push_line(out, &mut cur, q);
                }
                Ok(())
            }
            Kind::Line => {
                let tgt = cv.tgt.ok_or(Fail::Ring)?;
                let ps = tgt.param(s.t.unwrap_or(fv(s.p)));
                let pe = tgt.param(e.t.unwrap_or(fv(e.p)));
                let dir = (tgt.param(fv(run[m])) - tgt.param(fv(run[0]))).signum();
                if dir == 0.0 || dir * (pe - ps) <= 0.0 {
                    return Err(blame(&s, &e));
                }
                for &q in &run[1..m] {
                    let pq = tgt.param(fv(q));
                    let keep = (s.t.is_none() || dir * (pq - ps) > 1e-9)
                        && (e.t.is_none() || dir * (pe - pq) > 1e-9);
                    if keep {
                        push_line(out, &mut cur, q);
                    } else if hard(q) {
                        return Err(blame(&s, &e));
                    }
                }
                push_line(out, &mut cur, e.p);
                Ok(())
            }
            Kind::Circle(ci) => {
                let circ = &self.circles[ci];
                let tgt = cv.tgt.ok_or(Fail::Ring)?;
                let c = circ.c;
                // Unwrapped angles along the run.
                let mut th: Vec<f64> = Vec::with_capacity(run.len());
                th.push(ang(vsub(fv(run[0]), c)));
                for k in 1..=m {
                    let a = ang(vsub(fv(run[k]), c));
                    let prev = th[k - 1];
                    th.push(prev + norm_pi(a - prev));
                }
                let total = th[m] - th[0];
                let plain = s.t.is_none() && e.t.is_none();
                if total == 0.0 {
                    if plain {
                        for &q in &run[1..] {
                            push_line(out, &mut cur, q);
                        }
                        return Ok(());
                    }
                    return Err(blame(&s, &e));
                }
                let dir = total.signum();
                let near = |p: V, to: f64| -> f64 {
                    let a = ang(vsub(p, c));
                    to + norm_pi(a - to)
                };
                let ts = s.t.map_or(th[0], |t| near(t, th[0]));
                let te = e.t.map_or(th[m], |t| near(t, th[m]));
                if dir * (te - ts) <= 0.0 {
                    if plain {
                        for &q in &run[1..] {
                            push_line(out, &mut cur, q);
                        }
                        return Ok(());
                    }
                    return Err(blame(&s, &e));
                }
                // (point, angle, on the circle, hard)
                let mut ptv: Vec<(Point, f64, bool, bool)> = vec![(s.p, ts, s.t.is_some(), false)];
                for k in 1..m {
                    let keep = (s.t.is_none() || dir * (th[k] - ts) > 1e-12)
                        && (e.t.is_none() || dir * (te - th[k]) > 1e-12);
                    if keep {
                        ptv.push((run[k], th[k], false, hard(run[k])));
                    } else if hard(run[k]) {
                        return Err(blame(&s, &e));
                    }
                }
                ptv.push((e.p, te, e.t.is_some(), false));
                let exs = s.exempt.map(|(by, comp)| (by, shift_near(comp, ts)));
                let exe = e.exempt.map(|(by, comp)| (by, shift_near(comp, te)));
                let arcable = |i: usize| -> bool {
                    let (a, b) = (ptv[i], ptv[i + 1]);
                    !a.3 && !b.3
                        && dir * (b.1 - a.1) > 0.0
                        && clean(&cv.blocked, a.1.min(b.1), a.1.max(b.1), &tgt, exs, exe)
                };
                let proj = |th: f64| {
                    round_pt(c[0] + circ.r * libm::cos(th), c[1] + circ.r * libm::sin(th))
                };
                // A moved junction must be followed by its arc: a straight edge from the
                // true crossing to a vertex of `P` could cut through the curve.
                if s.moved && !arcable(0) {
                    return Err(Fail::Start);
                }
                if e.moved && !arcable(ptv.len() - 2) {
                    return Err(Fail::End);
                }
                let mut i = 0;
                while i + 1 < ptv.len() {
                    if arcable(i) {
                        let mut j = i + 1;
                        while j + 1 < ptv.len() && arcable(j) {
                            j += 1;
                        }
                        let (pa, ta, ona, _) = ptv[i];
                        let (pb, tb, onb, _) = ptv[j];
                        let a = if ona {
                            pa
                        } else {
                            proj(ta).map_err(|_| Fail::Ring)?
                        };
                        let b = if onb {
                            pb
                        } else {
                            proj(tb).map_err(|_| Fail::Ring)?
                        };
                        push_line(out, &mut cur, a);
                        self.emit_arc(circ, a, ta, b, tb, dir, out, &mut cur)
                            .ok_or(Fail::Ring)?;
                        push_line(out, &mut cur, pb);
                        i = j;
                    } else {
                        push_line(out, &mut cur, ptv[i + 1].0);
                        i += 1;
                    }
                }
                Ok(())
            }
        }
    }

    /// Appends arcs on `circ` from `a` (angle `ta`) to `b` (angle `tb`), split so no piece
    /// sweeps more than half a turn. Pieces too flat to be arcs become lines.
    #[allow(clippy::too_many_arguments)]
    fn emit_arc(
        &self,
        circ: &Circ,
        a: Point,
        ta: f64,
        b: Point,
        tb: f64,
        dir: f64,
        out: &mut Contour,
        cur: &mut Point,
    ) -> Option<()> {
        let sweep = dir * (tb - ta);
        if !sweep.is_finite() || sweep > 2.0 * TAU {
            return None;
        }
        let pieces = libm::ceil(sweep / PI - 1e-9).max(1.0) as usize;
        let mut prev = a;
        for k in 1..=pieces {
            let q = if k == pieces {
                b
            } else {
                let t = ta + (tb - ta) * k as f64 / pieces as f64;
                round_pt(
                    circ.c[0] + circ.r * libm::cos(t),
                    circ.c[1] + circ.r * libm::sin(t),
                )
                .ok()?
            };
            if q == prev {
                continue;
            }
            let s = sweep / pieces as f64;
            let sag = circ.r * (1.0 - libm::cos(s / 2.0));
            let co = fv(circ.c_out);
            let rp = vlen(vsub(fv(prev), co));
            // Direction implied by the rounded end points.
            let implied = {
                let d = ang(vsub(fv(q), co)) - ang(vsub(fv(prev), co));
                if dir > 0.0 { norm_tau(d) } else { norm_tau(-d) }
            };
            if sag < 0.5 || rp < 1.0 || prev == circ.c_out || q == circ.c_out {
                out.push(Curve::Line(q));
            } else if (implied - s).abs() > PI / 4.0 {
                return None;
            } else {
                out.push(Curve::CenterArc {
                    center: circ.c_out,
                    end: q,
                    ccw: dir > 0.0,
                });
            }
            *cur = q;
            prev = q;
        }
        Some(())
    }

    /// Whether a rebuilt contour is valid by construction: a single full circle, which is
    /// only produced when no other band comes near the circle.
    fn certain(&self, c: &Contour) -> bool {
        match c.as_slice() {
            [Curve::CenterArc { center, end, .. }] => dist2(*center, *end) >= 64,
            _ => false,
        }
    }

    /// Checks that a rebuilt ring's area is within the band bound of the polygon ring.
    fn area_ok(&self, c: &Contour, tr: &TaggedRing) -> bool {
        let n = tr.points.len();
        let mut per = 0.0;
        for i in 0..n {
            per += vlen(vsub(fv(tr.points[(i + 1) % n]), fv(tr.points[i])));
        }
        let ap = ring_area2(&tr.points) as f64 / 2.0;
        let ac = contour_area(c);
        (ac - ap).abs() <= per * (self.tol + 4.0) + 16.0
    }
}

/// A contour of straight lines through the points of a ring.
fn lines(pts: &[Point]) -> Contour {
    let n = pts.len();
    (0..n).map(|i| Curve::Line(pts[(i + 1) % n])).collect()
}

/// Rotates a contour so it starts at its lexicographically smallest vertex.
fn canonical_start(mut c: Contour) -> Contour {
    if let Some((k, _)) = c.iter().enumerate().min_by_key(|(i, e)| (e.end(), *i)) {
        let len = c.len();
        c.rotate_left((k + 1) % len);
    }
    c
}

/// Signed area of a curved contour (arcs as in [`Shape::to_polygon`]).
fn contour_area(c: &[Curve]) -> f64 {
    let Some(last) = c.last() else { return 0.0 };
    let o = last.end();
    let mut cur = o;
    let mut poly: i128 = 0;
    let mut extra = 0.0;
    for e in c {
        let end = e.end();
        poly += orient(o, cur, end);
        if let Curve::CenterArc { center, end, ccw } = *e
            && center != cur
        {
            let cc = fv(center);
            let r = vlen(vsub(fv(cur), cc));
            let a0 = ang(vsub(fv(cur), cc));
            let sweep = if end == cur {
                if ccw { TAU } else { -TAU }
            } else {
                sweep_to(a0, ang(vsub(fv(end), cc)), ccw)
            };
            extra += r * r * (sweep - libm::sin(sweep));
        }
        cur = end;
    }
    (poly as f64 + extra) / 2.0
}

/// Polygons (indices into `polys`) involved in a validity error.
fn error_polygons(e: &ValidityError, polys: &[Polygon]) -> Vec<usize> {
    let pts: Vec<V> = match e {
        ValidityError::TooFewVertices(r)
        | ValidityError::ZeroArea(r)
        | ValidityError::WrongOrientation(r) => return vec![r.polygon],
        ValidityError::DuplicateVertex(r, _) | ValidityError::SelfTouch(r, _) => {
            return vec![r.polygon];
        }
        ValidityError::SelfIntersection(p) => vec![[p.x, p.y]],
        ValidityError::OverlappingEdges(a, b) => vec![fv(*a), fv(*b)],
        ValidityError::InvalidNesting(p) | ValidityError::DisconnectedInterior(p) => vec![fv(*p)],
        ValidityError::CoordinateOutOfRange(p) => vec![fv(*p)],
        _ => return (0..polys.len()).collect(),
    };
    let mut hits = Vec::new();
    for (i, pg) in polys.iter().enumerate() {
        let near = pg.rings().any(|r| {
            let n = r.len();
            (0..n).any(|k| {
                let (a, b) = (fv(r[k]), fv(r[(k + 1) % n]));
                pts.iter().any(|&q| seg_dist(q, a, b) <= 2.0)
            })
        });
        if near {
            hits.push(i);
        }
    }
    if hits.is_empty() {
        (0..polys.len()).collect()
    } else {
        hits
    }
}

// ---------------------------------------------------------------------------------------
// Band separation geometry.

fn prim_bbox(p: &Prim) -> [f64; 4] {
    match *p {
        Prim::Seg { a, b } => [
            a[0].min(b[0]),
            a[1].min(b[1]),
            a[0].max(b[0]),
            a[1].max(b[1]),
        ],
        Prim::Arc {
            c,
            r,
            a0,
            sweep,
            e0,
            e1,
        } => {
            let mut bb = [
                e0[0].min(e1[0]),
                e0[1].min(e1[1]),
                e0[0].max(e1[0]),
                e0[1].max(e1[1]),
            ];
            for k in 0..4 {
                let t = k as f64 * PI / 2.0;
                if in_span(t, a0, sweep) {
                    let q = vadd(c, vmul(unit(t), r));
                    bb = [
                        bb[0].min(q[0]),
                        bb[1].min(q[1]),
                        bb[2].max(q[0]),
                        bb[3].max(q[1]),
                    ];
                }
            }
            bb
        }
    }
}

/// Whether angle `t` lies within the arc starting at `a0` with signed `sweep`.
fn in_span(t: f64, a0: f64, sweep: f64) -> bool {
    if sweep.abs() >= TAU - 1e-12 {
        return true;
    }
    if sweep >= 0.0 {
        norm_tau(t - a0) <= sweep
    } else {
        norm_tau(a0 - t) <= -sweep
    }
}

fn seg_dist(p: V, a: V, b: V) -> f64 {
    let d = vsub(b, a);
    let l2 = vdot(d, d);
    let t = if l2 > 0.0 {
        (vdot(vsub(p, a), d) / l2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    vlen(vsub(p, vadd(a, vmul(d, t))))
}

fn prim_dist(p: V, prim: &Prim) -> f64 {
    match *prim {
        Prim::Seg { a, b } => seg_dist(p, a, b),
        Prim::Arc {
            c,
            r,
            a0,
            sweep,
            e0,
            e1,
        } => {
            let q = vsub(p, c);
            let d = vlen(q);
            let ends = vlen(vsub(p, e0)).min(vlen(vsub(p, e1)));
            if d > 0.0 && in_span(ang(q), a0, sweep) {
                (d - r).abs().min(ends)
            } else if d == 0.0 {
                r.min(ends)
            } else {
                ends
            }
        }
    }
}

/// Parameters where the target meets the line through `p` with unit direction `d`.
fn tgt_line(t: &Tgt, p: V, d: V, out: &mut Vec<f64>) {
    match *t {
        Tgt::Circle { c, r } => {
            let w = vsub(c, p);
            let along = vdot(w, d);
            let h = vcross(d, w);
            let disc = (r - h.abs()) * (r + h.abs());
            if disc < 0.0 {
                return;
            }
            let s = libm::sqrt(disc);
            let foot = vsub(vadd(p, vmul(d, along)), c);
            for k in [-1.0, 1.0] {
                out.push(norm_tau(ang(vadd(foot, vmul(d, k * s)))));
            }
        }
        Tgt::Seg { a, d: td, .. } => {
            let den = vcross(td, d);
            if den.abs() < 1e-12 {
                return;
            }
            // a + s td on the line: cross(d, a + s td - p) = 0.
            let s = vcross(d, vsub(p, a)) / vcross(d, td);
            out.push(s);
        }
    }
}

/// Parameters where the target meets the circle `(q, rho)`.
fn tgt_circle(t: &Tgt, q: V, rho: f64, out: &mut Vec<f64>) {
    if rho.is_nan() || rho <= 0.0 {
        return;
    }
    match *t {
        Tgt::Circle { c, r } => {
            let dv = vsub(q, c);
            let dd = vlen(dv);
            if dd == 0.0 {
                return;
            }
            // Distance along the centre line from c: (dd^2 + r^2 - rho^2) / (2 dd).
            let a = (dd * dd + (r - rho) * (r + rho)) / (2.0 * dd);
            let h2 = (r - a) * (r + a);
            if h2 < 0.0 {
                return;
            }
            let h = libm::sqrt(h2);
            let base = ang(dv);
            for k in [-1.0, 1.0] {
                out.push(norm_tau(base + libm::atan2(k * h, a)));
            }
        }
        Tgt::Seg { a, d, .. } => {
            let w = vsub(q, a);
            let along = vdot(w, d);
            let h = vcross(d, w);
            let disc = (rho - h.abs()) * (rho + h.abs());
            if disc < 0.0 {
                return;
            }
            let s = libm::sqrt(disc);
            out.push(along - s);
            out.push(along + s);
        }
    }
}

/// Parameter intervals of the target lying within `wsum` of the primitive.
fn blocked_intervals(t: &Tgt, prim: &Prim, wsum: f64) -> Vec<(f64, f64)> {
    let mut cands: Vec<f64> = Vec::new();
    match *prim {
        Prim::Seg { a, b } => {
            let len = vlen(vsub(b, a));
            if len > 0.0 {
                let d = vmul(vsub(b, a), 1.0 / len);
                let nrm = [-d[1], d[0]];
                tgt_line(t, vadd(a, vmul(nrm, wsum)), d, &mut cands);
                tgt_line(t, vsub(a, vmul(nrm, wsum)), d, &mut cands);
            }
            tgt_circle(t, a, wsum, &mut cands);
            tgt_circle(t, b, wsum, &mut cands);
        }
        Prim::Arc { c, r, e0, e1, .. } => {
            tgt_circle(t, c, r + wsum, &mut cands);
            tgt_circle(t, c, r - wsum, &mut cands);
            tgt_circle(t, e0, wsum, &mut cands);
            tgt_circle(t, e1, wsum, &mut cands);
        }
    }
    let (hi, periodic) = match *t {
        Tgt::Circle { .. } => (TAU, true),
        Tgt::Seg { len, .. } => (len, false),
    };
    cands.retain(|x| x.is_finite() && *x > 0.0 && *x < hi);
    cands.push(0.0);
    if !periodic {
        cands.push(hi);
    }
    cands.sort_by(f64::total_cmp);
    cands.dedup();
    let mut spans: Vec<(f64, f64)> = cands.windows(2).map(|w| (w[0], w[1])).collect();
    if periodic {
        let first = cands[0];
        let last = *cands.last().unwrap_or(&0.0);
        spans.push((last, first + TAU));
    }
    let eps = t.eps();
    let mut res: Vec<(f64, f64)> = Vec::new();
    for (lo, h) in spans {
        let mid = (lo + h) / 2.0;
        if prim_dist(t.at(mid), prim) <= wsum {
            match res.last_mut() {
                Some(l) if l.1 >= lo - 1e-12 => l.1 = h,
                _ => res.push((lo, h)),
            }
        }
    }
    let mut out = Vec::new();
    for (lo, h) in res {
        let (lo, h) = (lo - eps, h + eps);
        if periodic {
            if h - lo >= TAU {
                out.push((0.0, TAU));
            } else if lo < 0.0 {
                out.push((lo + TAU, TAU));
                out.push((0.0, h));
            } else if h > TAU {
                out.push((lo, TAU));
                out.push((0.0, h - TAU));
            } else {
                out.push((lo, h));
            }
        } else {
            out.push((lo.max(0.0), h.min(hi)));
        }
    }
    out
}

/// Shifts a parameter of a periodic target by whole turns to lie closest to `comp`.
fn shift_into(p: f64, comp: (f64, f64), t: &Tgt) -> f64 {
    if !t.periodic() {
        return p;
    }
    let mid = (comp.0 + comp.1) / 2.0;
    p + TAU * libm::round((mid - p) / TAU)
}

/// Shifts an interval by whole turns so its middle is closest to `to`.
fn shift_near(comp: (f64, f64), to: f64) -> (f64, f64) {
    let mid = (comp.0 + comp.1) / 2.0;
    let k = TAU * libm::round((to - mid) / TAU);
    (comp.0 + k, comp.1 + k)
}

fn inside(p: f64, comp: (f64, f64), eps: f64) -> bool {
    p >= comp.0 - eps && p <= comp.1 + eps
}

fn shifts(t: &Tgt) -> &'static [f64] {
    if t.periodic() {
        &[-2.0, -1.0, 0.0, 1.0, 2.0]
    } else {
        &[0.0]
    }
}

/// Connected union of the intervals caused by `by` that contains parameter `q`. `None`
/// when there is none or it covers a whole turn.
fn component(list: &[(f64, f64, u32)], by: u32, q: f64, t: &Tgt) -> Option<(f64, f64)> {
    let eps = t.eps();
    let mut comp: Option<(f64, f64)> = None;
    loop {
        let mut changed = false;
        for &(lo, hi, b) in list {
            if b != by {
                continue;
            }
            for &k in shifts(t) {
                let (lo, hi) = (lo + k * TAU, hi + k * TAU);
                match comp {
                    None => {
                        if q >= lo - eps && q <= hi + eps {
                            comp = Some((lo, hi));
                            changed = true;
                        }
                    }
                    Some((cl, ch)) => {
                        if lo <= ch + 1e-12 && hi >= cl - 1e-12 && (lo < cl || hi > ch) {
                            comp = Some((cl.min(lo), ch.max(hi)));
                            changed = true;
                        }
                    }
                }
            }
        }
        let c = comp?;
        if t.periodic() && c.1 - c.0 >= TAU - 1e-9 {
            return None;
        }
        if !changed {
            return comp;
        }
    }
}

/// Whether an interval caused by a curve other than `by` overlaps `comp`.
fn others_overlap(list: &[(f64, f64, u32)], by: u32, comp: (f64, f64), t: &Tgt) -> bool {
    list.iter().any(|&(lo, hi, b)| {
        b != by
            && shifts(t)
                .iter()
                .any(|&k| lo + k * TAU <= comp.1 + 1e-12 && hi + k * TAU >= comp.0 - 1e-12)
    })
}

/// Whether `[lo, hi]` (unwrapped angles of a circle, or line parameters) is free of
/// blocked intervals, ignoring those of the resolved junction components.
fn clean(
    list: &[(f64, f64, u32)],
    lo: f64,
    hi: f64,
    t: &Tgt,
    ex_a: Option<(u32, (f64, f64))>,
    ex_b: Option<(u32, (f64, f64))>,
) -> bool {
    if !lo.is_finite() || !hi.is_finite() || hi - lo > 4.0 * TAU {
        return false;
    }
    let kmin = libm::floor((lo - TAU) / TAU) as i64;
    let kmax = libm::ceil(hi / TAU) as i64;
    for &(il, ih, b) in list {
        let ks: Vec<f64> = if t.periodic() {
            (kmin..=kmax).map(|k| k as f64).collect()
        } else {
            vec![0.0]
        };
        for k in ks {
            let (l, h) = (il + k * TAU, ih + k * TAU);
            if l > hi + 1e-12 || h < lo - 1e-12 {
                continue;
            }
            let exempt = [ex_a, ex_b]
                .iter()
                .flatten()
                .any(|&(by, comp)| by == b && l <= comp.1 + 1e-12 && h >= comp.0 - 1e-12);
            if !exempt {
                return false;
            }
        }
    }
    true
}

/// Intersections of two curves (circles in full, segments as infinite lines), with the
/// sine of the crossing angle.
fn intersections(a: &Tgt, b: &Tgt) -> Vec<(V, f64)> {
    let mut out = Vec::new();
    match (*a, *b) {
        (Tgt::Circle { c: c1, r: r1 }, Tgt::Circle { c: c2, r: r2 }) => {
            let mut ts = Vec::new();
            tgt_circle(a, c2, r2, &mut ts);
            for t in ts {
                let p = vadd(c1, vmul(unit(t), r1));
                let sin = (vcross(vsub(p, c1), vsub(p, c2)) / (r1 * r2)).abs();
                out.push((p, sin));
            }
        }
        (Tgt::Circle { c, r }, Tgt::Seg { a: la, d, .. })
        | (Tgt::Seg { a: la, d, .. }, Tgt::Circle { c, r }) => {
            let mut ts = Vec::new();
            tgt_line(&Tgt::Circle { c, r }, la, d, &mut ts);
            for t in ts {
                let p = vadd(c, vmul(unit(t), r));
                let sin = (vdot(d, vsub(p, c)) / r).abs();
                out.push((p, sin));
            }
        }
        _ => {}
    }
    out
}

#[cfg(test)]
#[allow(clippy::cloned_ref_to_slice_refs)]
mod tests {
    use super::*;
    use crate::geom::Ring;

    fn p(x: i64, y: i64) -> Point {
        Point::new(x, y)
    }

    fn circle(c: Point, r: i64) -> Shape {
        Shape::new(
            vec![Curve::CenterArc {
                center: c,
                end: p(c.x + r, c.y),
                ccw: true,
            }],
            vec![],
        )
    }

    fn rounded_rect(x0: i64, y0: i64, x1: i64, y1: i64, r: i64) -> Shape {
        Shape::new(
            vec![
                Curve::Line(p(x1 - r, y0)),
                Curve::CenterArc {
                    center: p(x1 - r, y0 + r),
                    end: p(x1, y0 + r),
                    ccw: true,
                },
                Curve::Line(p(x1, y1 - r)),
                Curve::CenterArc {
                    center: p(x1 - r, y1 - r),
                    end: p(x1 - r, y1),
                    ccw: true,
                },
                Curve::Line(p(x0 + r, y1)),
                Curve::CenterArc {
                    center: p(x0 + r, y1 - r),
                    end: p(x0, y1 - r),
                    ccw: true,
                },
                Curve::Line(p(x0, y0 + r)),
                Curve::CenterArc {
                    center: p(x0 + r, y0 + r),
                    end: p(x0 + r, y0),
                    ccw: true,
                },
            ],
            vec![],
        )
    }

    fn rect(x0: i64, y0: i64, x1: i64, y1: i64) -> Shape {
        Shape::new(
            vec![
                Curve::Line(p(x1, y0)),
                Curve::Line(p(x1, y1)),
                Curve::Line(p(x0, y1)),
                Curve::Line(p(x0, y0)),
            ],
            vec![],
        )
    }

    fn arcs(s: &[Shape]) -> Vec<(Point, Point, Point, bool)> {
        let mut v = Vec::new();
        for sh in s {
            for c in core::iter::once(&sh.contour).chain(sh.holes.iter()) {
                let mut cur = c.last().map(|e| e.end()).unwrap_or_default();
                for e in c {
                    if let Curve::CenterArc { center, end, ccw } = *e {
                        v.push((cur, center, end, ccw));
                    }
                    cur = e.end();
                }
            }
        }
        v
    }

    fn approx(s: &[Shape], t: i64) -> Vec<Polygon> {
        s.iter()
            .map(|x| x.to_polygon(ArcTol::new(t, Side::Nearest)).unwrap())
            .collect()
    }

    fn area(s: &[Shape]) -> f64 {
        s.iter()
            .map(|x| {
                contour_area(&x.contour) + x.holes.iter().map(|h| contour_area(h)).sum::<f64>()
            })
            .sum()
    }

    fn check(s: &[Shape]) {
        validate_set(&approx(s, 1)).unwrap();
        validate_set(&approx(s, 50)).unwrap();
    }

    fn on_circle(out: &[Shape], c: Point, r: f64) {
        on_circles(out, c, &[r]);
    }

    fn on_circles(out: &[Shape], c: Point, radii: &[f64]) {
        for (st, _, en, _) in arcs(out).into_iter().filter(|a| a.1 == c) {
            let r = vlen(vsub(fv(st), fv(c)));
            assert!(radii.iter().any(|x| (x - r).abs() <= 0.75), "{st:?} at {r}");
            let d = vlen(vsub(fv(en), fv(c)));
            assert!((d - r).abs() <= 1.5, "{en:?} at {d}, expected {r}");
        }
    }

    const T: ArcTol = ArcTol::new(100, Side::Inside);

    #[test]
    fn circle_union_circle() {
        let a = circle(p(0, 0), 1_000_000);
        let b = circle(p(1_200_000, 0), 1_000_000);
        for side in [Side::Inside, Side::Outside, Side::Nearest] {
            let out = curved_boolean(
                Op::Union,
                &[a.clone()],
                &[b.clone()],
                FillRule::NonZero,
                ArcTol::new(100, side),
            )
            .unwrap();
            assert_eq!(out.len(), 1);
            assert!(out[0].holes.is_empty());
            // Two arcs (one per circle), meeting at the rounded true crossings.
            let a = arcs(&out);
            assert!(
                out[0]
                    .contour
                    .iter()
                    .all(|e| matches!(e, Curve::CenterArc { .. }))
            );
            assert_eq!(a.len(), 4, "{:?}", out[0].contour); // each arc split in halves
            let ys: Vec<i64> = out[0].contour.iter().map(|e| e.end().y.abs()).collect();
            // True crossings at x = 600 000, y = +-800 000.
            assert!(ys.contains(&800_000), "{ys:?}");
            on_circle(&out, p(0, 0), 1e6);
            on_circle(&out, p(1_200_000, 0), 1e6);
            check(&out);
            // Exact area of the union of two unit disks at distance 1.2.
            let r: f64 = 1e6;
            let d: f64 = 1.2e6;
            let lens =
                2.0 * r * r * libm::acos(d / (2.0 * r)) - d / 2.0 * libm::sqrt(4.0 * r * r - d * d);
            let exact = 2.0 * PI * r * r - lens;
            assert!(
                (area(&out) - exact).abs() < 1e-6 * exact,
                "{} {exact}",
                area(&out)
            );
        }
    }

    #[test]
    fn rounded_rect_minus_circle() {
        let mm = 1_000_000;
        let board = rounded_rect(0, 0, 10 * mm, 8 * mm, mm);
        // One via inside, one crossing the right edge, one crossing the corner arc.
        let vias = vec![
            circle(p(5 * mm, 4 * mm), 300_000),
            circle(p(10 * mm, 4 * mm), 500_000),
            circle(p(9_700_000, 7_700_000), 400_000),
        ];
        let out = curved_boolean(Op::Difference, &[board], &vias, FillRule::NonZero, T).unwrap();
        check(&out);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].holes.len(), 1);
        assert_eq!(
            out[0].holes[0],
            vec![Curve::CenterArc {
                center: p(5 * mm, 4 * mm),
                end: p(5_300_000, 4 * mm),
                ccw: false,
            }]
        );
        let a = arcs(&out);
        // Three untouched corner arcs (each a quarter), the notch arcs.
        for c in [p(mm, mm), p(9 * mm, mm), p(mm, 7 * mm)] {
            assert_eq!(
                a.iter().filter(|x| x.1 == c).count(),
                1,
                "{:?}",
                out[0].contour
            );
        }
        assert!(a.iter().any(|x| x.1 == p(10 * mm, 4 * mm)));
        assert!(a.iter().any(|x| x.1 == p(9_700_000, 7_700_000)));
        assert!(a.iter().any(|x| x.1 == p(9 * mm, 7 * mm)));
        // Notch crossing the right edge: ends at x = 10 mm, y = 4 mm +- 0.5 mm.
        assert!(
            out[0]
                .contour
                .iter()
                .any(|e| e.end() == p(10 * mm, 4_500_000))
        );
        assert!(
            out[0]
                .contour
                .iter()
                .any(|e| e.end() == p(10 * mm, 3_500_000))
        );
        on_circle(&out, p(9 * mm, 7 * mm), mm as f64);
        on_circle(&out, p(9_700_000, 7_700_000), 400_000.0);
        // Few elements: 4 lines + 4 corner arcs (+ splits) + 2 notches.
        assert!(out[0].contour.len() <= 16, "{:?}", out[0].contour);
    }

    #[test]
    fn annulus_and_concentric() {
        let outer = circle(p(0, 0), 1_000_000);
        let inner = circle(p(0, 0), 400_000);
        let out = curved_boolean(
            Op::Difference,
            &[outer.clone()],
            &[inner.clone()],
            FillRule::NonZero,
            T,
        )
        .unwrap();
        check(&out);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].contour.len(), 1);
        assert_eq!(out[0].holes.len(), 1);
        assert_eq!(out[0].holes[0].len(), 1);
        let exact = PI * (1e12 - 0.16e12);
        assert!((area(&out) - exact).abs() < 1e-9 * exact);
        // Annulus shape as input (hole) intersected with a concentric disk.
        let ann = Shape::new(outer.contour.clone(), vec![inner.contour.clone()]);
        let mid = circle(p(0, 0), 700_000);
        let out = curved_boolean(Op::Intersection, &[ann], &[mid], FillRule::NonZero, T).unwrap();
        check(&out);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].contour.len(), 1);
        assert_eq!(out[0].holes[0].len(), 1);
        on_circles(&out, p(0, 0), &[400_000.0, 700_000.0]);
        // Union of concentric circles is the big one.
        let out =
            curved_boolean(Op::Union, &[outer.clone()], &[inner], FillRule::NonZero, T).unwrap();
        assert_eq!(out, vec![outer_canon(&outer)]);
    }

    fn outer_canon(s: &Shape) -> Shape {
        Shape::new(canonical_start(s.contour.clone()), vec![])
    }

    #[test]
    fn tangent_circles() {
        // Externally tangent: the union keeps both circles as arcs away from the contact,
        // and stays valid.
        let a = circle(p(0, 0), 1_000_000);
        let b = circle(p(2_000_000, 0), 1_000_000);
        for side in [Side::Inside, Side::Outside, Side::Nearest] {
            let tol = ArcTol::new(100, side);
            let out = curved_boolean(
                Op::Union,
                &[a.clone()],
                &[b.clone()],
                FillRule::NonZero,
                tol,
            )
            .unwrap();
            check(&out);
            let a = arcs(&out);
            assert!(a.iter().any(|x| x.1 == p(0, 0)));
            assert!(a.iter().any(|x| x.1 == p(2_000_000, 0)));
            let exact = 2.0 * PI * 1e12;
            assert!((area(&out) - exact).abs() < 1e-3 * exact);
            // Internally tangent difference.
            let c = circle(p(500_000, 0), 500_000);
            let out = curved_boolean(
                Op::Difference,
                &[circle(p(0, 0), 1_000_000)],
                &[c],
                FillRule::NonZero,
                tol,
            )
            .unwrap();
            check(&out);
            let exact = PI * (1e12 - 0.25e12);
            assert!((area(&out) - exact).abs() < 1e-3 * exact);
        }
    }

    #[test]
    fn coincident_arcs() {
        // Same circle written twice with different start points: one circle out.
        let a = circle(p(0, 0), 1_000_000);
        let b = Shape::new(
            vec![
                Curve::CenterArc {
                    center: p(0, 0),
                    end: p(-1_000_000, 0),
                    ccw: true,
                },
                Curve::CenterArc {
                    center: p(0, 0),
                    end: p(1_000_000, 0),
                    ccw: true,
                },
            ],
            vec![],
        );
        for side in [Side::Inside, Side::Outside, Side::Nearest] {
            let tol = ArcTol::new(100, side);
            let u = curved_boolean(
                Op::Union,
                &[a.clone()],
                &[b.clone()],
                FillRule::NonZero,
                tol,
            )
            .unwrap();
            assert_eq!(u.len(), 1);
            assert_eq!(u[0].contour.len(), 1, "{side:?} {:?}", u[0].contour);
            let i = curved_boolean(
                Op::Intersection,
                &[a.clone()],
                &[b.clone()],
                FillRule::NonZero,
                tol,
            )
            .unwrap();
            assert_eq!(i.len(), 1);
            assert_eq!(i[0].contour.len(), 1);
            let d = curved_boolean(
                Op::Difference,
                &[a.clone()],
                &[b.clone()],
                FillRule::NonZero,
                tol,
            )
            .unwrap();
            assert!(d.is_empty(), "{side:?} {d:?}");
            let x = curved_boolean(Op::Xor, &[a.clone()], &[b.clone()], FillRule::NonZero, tol)
                .unwrap();
            assert!(x.is_empty(), "{side:?} {x:?}");
        }
        // Two rounded rects sharing a corner arc: union is one shape with arcs.
        let mm = 1_000_000;
        let r1 = rounded_rect(0, 0, 4 * mm, 4 * mm, mm);
        let r2 = rounded_rect(0, 0, 4 * mm, 6 * mm, mm);
        let u = curved_boolean(Op::Union, &[r1], &[r2], FillRule::NonZero, T).unwrap();
        check(&u);
        assert_eq!(u.len(), 1);
        assert_eq!(arcs(&u).len(), 4);
    }

    #[test]
    fn grazing_lines() {
        // A long thin wedge crossing a circle at a grazing angle.
        let c = circle(p(0, 0), 1_000_000);
        for dy in [999_000, 999_990, 1_000_000, 1_000_010, 1_001_000] {
            let wedge = Shape::new(
                vec![
                    Curve::Line(p(3_000_000, dy)),
                    Curve::Line(p(3_000_000, dy + 1_000)),
                    Curve::Line(p(-3_000_000, dy + 3_000)),
                    Curve::Line(p(-3_000_000, dy)),
                ],
                vec![],
            );
            for op in [Op::Union, Op::Difference, Op::Intersection, Op::Xor] {
                for side in [Side::Inside, Side::Outside, Side::Nearest] {
                    let out = curved_boolean(
                        op,
                        &[c.clone()],
                        &[wedge.clone()],
                        FillRule::NonZero,
                        ArcTol::new(100, side),
                    )
                    .unwrap();
                    check(&out);
                    on_circle(&out, p(0, 0), 1e6);
                }
            }
        }
    }

    #[test]
    fn steep_line_crossing_is_exact() {
        // A rectangle edge crossing a circle at a 45 degree angle: the junction is the
        // rounded true crossing.
        let c = circle(p(0, 0), 1_000_000);
        let r = rect(0, 0, 3_000_000, 3_000_000);
        let out = curved_boolean(Op::Difference, &[c], &[r], FillRule::NonZero, T).unwrap();
        check(&out);
        let ends: Vec<Point> = out[0].contour.iter().map(|e| e.end()).collect();
        assert!(ends.contains(&p(1_000_000, 0)), "{ends:?}");
        assert!(ends.contains(&p(0, 1_000_000)), "{ends:?}");
        assert!(ends.contains(&p(0, 0)), "{ends:?}");
        let exact = 0.75 * PI * 1e12;
        assert!(
            (area(&out) - exact).abs() < 1e-5 * exact,
            "{} {exact}",
            area(&out)
        );
    }

    #[test]
    fn tiny_and_huge_radii() {
        for r in [1, 2, 3, 5, 10] {
            let a = circle(p(0, 0), r);
            let b = circle(p(r, 0), r);
            for op in [Op::Union, Op::Difference, Op::Intersection, Op::Xor] {
                let out = curved_boolean(
                    op,
                    &[a.clone()],
                    &[b.clone()],
                    FillRule::NonZero,
                    ArcTol::new(1, Side::Nearest),
                )
                .unwrap();
                validate_set(&approx(&out, 1)).unwrap();
            }
        }
        // Huge circles near the coordinate limit.
        let big = 1i64 << 39;
        let a = circle(p(-(1 << 38), 0), big);
        let b = circle(p(1 << 38, 0), big);
        for side in [Side::Inside, Side::Nearest] {
            let tol = ArcTol::new(1 << 20, side);
            let out = curved_boolean(
                Op::Intersection,
                &[a.clone()],
                &[b.clone()],
                FillRule::NonZero,
                tol,
            )
            .unwrap();
            assert_eq!(out.len(), 1);
            validate_set(&approx(&out, 1 << 20)).unwrap();
            on_circle(&out, p(-(1 << 38), 0), big as f64);
            on_circle(&out, p(1 << 38, 0), big as f64);
            assert!(
                out[0]
                    .contour
                    .iter()
                    .all(|e| matches!(e, Curve::CenterArc { .. }))
            );
        }
        // A huge arc (outline of a huge rounded rect) minus a small circle on it.
        let m = (1i64 << 40) - 2_000_000;
        let s = rounded_rect(-m, -m, m, m, 1 << 39);
        let hole = circle(p(m, 0), 1_000_000);
        let out = curved_boolean(
            Op::Difference,
            &[s],
            &[hole],
            FillRule::NonZero,
            ArcTol::new(1_000, Side::Inside),
        )
        .unwrap();
        assert_eq!(out.len(), 1);
        assert_eq!(arcs(&out).iter().filter(|x| x.1 == p(m, 0)).count(), 1);
    }

    #[test]
    fn three_point_arcs_and_degenerate_input() {
        // Half disk via a three-point arc with a non-integral centre.
        let s = Shape::new(
            vec![
                Curve::Arc {
                    mid: p(500_001, 500_000),
                    end: p(1_000_001, 0),
                },
                Curve::Line(p(0, 0)),
            ],
            vec![],
        );
        let r = rect(400_000, -100_000, 600_000, 700_000);
        let out = curved_boolean(Op::Difference, &[s], &[r], FillRule::NonZero, T).unwrap();
        check(&out);
        assert_eq!(out.len(), 2);
        assert!(!arcs(&out).is_empty());
        // Empty and degenerate input.
        assert!(
            curved_boolean(Op::Union, &[], &[], FillRule::NonZero, T)
                .unwrap()
                .is_empty()
        );
        let deg = Shape::new(vec![Curve::Line(p(0, 0))], vec![]);
        assert!(
            curved_boolean(Op::Union, &[deg], &[], FillRule::NonZero, T)
                .unwrap()
                .is_empty()
        );
        assert!(
            curved_boolean(
                Op::Union,
                &[],
                &[],
                FillRule::NonZero,
                ArcTol::new(0, Side::Inside)
            )
            .is_err()
        );
        let far = circle(p(1 << 41, 0), 10);
        assert!(curved_boolean(Op::Union, &[far], &[], FillRule::NonZero, T).is_err());
        // Off-circle centre-form arc: kept as its approximation, still valid.
        let off = Shape::new(
            vec![
                Curve::CenterArc {
                    center: p(0, 0),
                    end: p(0, 1_000_100),
                    ccw: true,
                },
                Curve::Line(p(1_000_000, 0)),
            ],
            vec![],
        );
        let out = curved_boolean(Op::Union, &[off], &[], FillRule::NonZero, T).unwrap();
        check(&out);
        assert!(arcs(&out).is_empty());
    }

    #[test]
    fn side_guarantee_zone_minus_obstacles() {
        // A circle grazing the zone's rounded corner from outside by less than the
        // tolerance: the Inside result must not intrude into it.
        let mm = 1_000_000;
        let zone = rounded_rect(0, 0, 10 * mm, 10 * mm, 2 * mm);
        // Corner centre (8, 8) mm, radius 2 mm; diagonal point at distance 2 mm.
        let d = (2.0 * mm as f64 + 400_000.0 - 50.0) / 2f64.sqrt();
        let c = p(8 * mm + d as i64, 8 * mm + d as i64);
        let obstacle = circle(c, 400_000);
        let tol = ArcTol::new(2_000, Side::Inside);
        let out =
            curved_boolean(Op::Difference, &[zone], &[obstacle], FillRule::NonZero, tol).unwrap();
        check(&out);
        let poly = approx(&out, 1);
        let fine = &poly[0];
        // Points just inside the obstacle are not in the result.
        for k in 0..64 {
            let a = k as f64 / 64.0 * TAU;
            let q = p(
                c.x + (libm::cos(a) * 399_990.0) as i64,
                c.y + (libm::sin(a) * 399_990.0) as i64,
            );
            assert_ne!(
                crate::query::locate_in_polygon(fine, q),
                crate::query::Location::Inside,
                "{q:?}"
            );
        }
        let _ = Ring::default();
    }
}
