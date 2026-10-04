//! Boolean operations on polygons (and clipping of open paths by polygons).

use crate::arrangement::{Arrangement, Engine, InEdge, Noding};
use crate::assemble::{DirEdge, assemble};
use crate::error::{Error, Result};
use crate::geom::{
    Path, Point, PolyTree, Polygon, PolygonSet, Ring, TaggedPath, TaggedPolygon, TaggedRing,
};

/// Rule deciding which regions are "inside" from their winding number.
///
/// The winding number of a point counts how many times the input rings wind around it,
/// counter-clockwise rings counting `+1` (Y-up convention).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum FillRule {
    /// Odd winding numbers are inside.
    EvenOdd,
    /// Non-zero winding numbers are inside (the usual rule for unions of shapes).
    #[default]
    NonZero,
    /// Strictly positive winding numbers are inside.
    Positive,
    /// Strictly negative winding numbers are inside.
    Negative,
}

impl FillRule {
    /// Whether a winding number is inside under this rule.
    #[inline]
    pub fn is_inside(self, w: i32) -> bool {
        match self {
            FillRule::EvenOdd => w & 1 != 0,
            FillRule::NonZero => w != 0,
            FillRule::Positive => w > 0,
            FillRule::Negative => w < 0,
        }
    }
}

/// Boolean operation.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Op {
    /// Subject ∪ clip.
    #[default]
    Union,
    /// Subject ∩ clip.
    Intersection,
    /// Subject − clip.
    Difference,
    /// Subject ⊕ clip (symmetric difference).
    Xor,
}

impl Op {
    #[inline]
    fn apply(self, s: bool, c: bool) -> bool {
        match self {
            Op::Union => s || c,
            Op::Intersection => s && c,
            Op::Difference => s && !c,
            Op::Xor => s != c,
        }
    }
}

/// Callback receiving a vertex sequence and its optional per-edge tags.
pub type VertexVisitor<'a> = dyn FnMut(&[Point], Option<&[u64]>) + 'a;

/// Anything that can be fed to a boolean operation as a set of closed rings.
///
/// Implemented for [`Ring`], [`Polygon`], [`TaggedRing`], [`TaggedPolygon`], [`PolyTree`],
/// raw vertex vectors, and slices / `Vec`s / arrays of any of those. Rings may have any
/// orientation, self-intersect and overlap; the fill rule decides what is inside.
pub trait RingSource {
    /// Calls `f(points, tags)` for every ring. `tags`, when present, has one tag per edge.
    fn visit_rings(&self, f: &mut VertexVisitor<'_>);

    /// When the source is made only of [`Polygon`]s (untagged), calls `f` on each and
    /// returns `true`; otherwise returns `false` (possibly after some calls). Lets operations
    /// detect already canonical input and skip re-normalizing it.
    fn visit_polygons<'a>(&'a self, _f: &mut dyn FnMut(&'a Polygon)) -> bool {
        false
    }
}

impl RingSource for Ring {
    fn visit_rings(&self, f: &mut VertexVisitor<'_>) {
        f(&self.0, None)
    }
}

impl RingSource for Vec<Point> {
    fn visit_rings(&self, f: &mut VertexVisitor<'_>) {
        f(self, None)
    }
}

impl RingSource for TaggedRing {
    fn visit_rings(&self, f: &mut VertexVisitor<'_>) {
        f(&self.points, Some(&self.tags))
    }
}

impl RingSource for Polygon {
    fn visit_rings(&self, f: &mut VertexVisitor<'_>) {
        for r in self.rings() {
            f(&r.0, None)
        }
    }
    fn visit_polygons<'a>(&'a self, f: &mut dyn FnMut(&'a Polygon)) -> bool {
        f(self);
        true
    }
}

impl RingSource for TaggedPolygon {
    fn visit_rings(&self, f: &mut VertexVisitor<'_>) {
        f(&self.outer.points, Some(&self.outer.tags));
        for h in &self.holes {
            f(&h.points, Some(&h.tags));
        }
    }
}

impl RingSource for PolyTree {
    fn visit_rings(&self, f: &mut VertexVisitor<'_>) {
        for n in &self.nodes {
            f(&n.ring.0, Some(&n.tags))
        }
    }
}

impl<T: RingSource> RingSource for [T] {
    fn visit_rings(&self, f: &mut VertexVisitor<'_>) {
        for x in self {
            x.visit_rings(f)
        }
    }
    fn visit_polygons<'a>(&'a self, f: &mut dyn FnMut(&'a Polygon)) -> bool {
        self.iter().all(|x| x.visit_polygons(f))
    }
}

impl<T: RingSource> RingSource for Vec<T> {
    fn visit_rings(&self, f: &mut VertexVisitor<'_>) {
        self.as_slice().visit_rings(f)
    }
    fn visit_polygons<'a>(&'a self, f: &mut dyn FnMut(&'a Polygon)) -> bool {
        self.as_slice().visit_polygons(f)
    }
}

impl<T: RingSource, const N: usize> RingSource for [T; N] {
    fn visit_rings(&self, f: &mut VertexVisitor<'_>) {
        self.as_slice().visit_rings(f)
    }
    fn visit_polygons<'a>(&'a self, f: &mut dyn FnMut(&'a Polygon)) -> bool {
        self.as_slice().visit_polygons(f)
    }
}

impl<T: RingSource + ?Sized> RingSource for &T {
    fn visit_rings(&self, f: &mut VertexVisitor<'_>) {
        (**self).visit_rings(f)
    }
    fn visit_polygons<'a>(&'a self, f: &mut dyn FnMut(&'a Polygon)) -> bool {
        (**self).visit_polygons(f)
    }
}

/// Reusable boolean-operation engine.
///
/// ```
/// use polyclip::{Boolean, FillRule, Op, Ring};
/// let a = Ring::from([(0, 0), (10, 0), (10, 10), (0, 10)]);
/// let b = Ring::from([(5, 5), (15, 5), (15, 15), (5, 15)]);
/// let u = Boolean::new()
///     .subject(&a, FillRule::NonZero)
///     .clip(&b, FillRule::NonZero)
///     .op(Op::Union)
///     .execute()
///     .unwrap();
/// assert_eq!(u.len(), 1);
/// assert_eq!(u[0].outer.len(), 8);
/// ```
///
/// # Semantics
///
/// * Inputs are rings in any orientation; each operand has its own [`FillRule`].
/// * The result is computed on the snap-rounded arrangement of all input edges: every
///   intersection vertex is rounded to the nearest integer point (`floor(v + 1/2)` per
///   coordinate), and all edges passing within the same unit pixel are routed through that
///   point. No point moves by more than `sqrt(2)/2` units, and the output is always valid:
///   simple rings, no crossing edges, correct nesting.
/// * Output is canonical: outer rings counter-clockwise, holes clockwise, every ring starting
///   at its lexicographically smallest vertex, no repeated vertices, no zero-length edges,
///   collinear vertices removed (unless [`keep_collinear`](Self::keep_collinear)), polygons
///   sorted by outer ring vertex sequence and holes by their vertex sequence.
/// * Regions touching at a single point are separate rings: two polygons meeting at a corner
///   are two polygons, and a hole touching its outer ring at a vertex is a hole sharing that
///   vertex. Rings never touch along an edge.
/// * Degenerate input (zero-length edges, rings with fewer than three vertices, zero-area
///   rings, duplicates) is accepted; it simply contributes nothing where it has no area.
/// * Edge tags: every output edge carries the tag of the input edge it lies on. When several
///   input edges coincide, the tag comes from the operand that matters there (non-zero net
///   winding change), preferring the subject, then the smallest tag. Untagged rings carry
///   tag `0`.
#[derive(Clone, Debug, Default)]
pub struct Boolean {
    edges: Vec<InEdge>,
    /// Index of the first edge of every ring in `edges` (rings are contiguous).
    ring_starts: Vec<u32>,
    fill: [FillRule; 2],
    op: Op,
    keep_collinear: bool,
    clustering: Clustering,
    error: Option<Error>,
}

impl Boolean {
    /// Creates an empty engine (union, non-zero fill rules).
    pub fn new() -> Self {
        Self::default()
    }

    /// Removes all input, keeping allocations and settings.
    pub fn clear(&mut self) {
        self.edges.clear();
        self.ring_starts.clear();
        self.error = None;
    }

    /// Adds subject rings and sets the subject fill rule (builder form).
    pub fn subject(mut self, rings: &(impl RingSource + ?Sized), rule: FillRule) -> Self {
        self.add_subject(rings, rule);
        self
    }

    /// Adds clip rings and sets the clip fill rule (builder form).
    pub fn clip(mut self, rings: &(impl RingSource + ?Sized), rule: FillRule) -> Self {
        self.add_clip(rings, rule);
        self
    }

    /// Sets the operation (builder form).
    pub fn op(mut self, op: Op) -> Self {
        self.op = op;
        self
    }

    /// Keeps collinear vertices in the output (builder form). Off by default.
    pub fn keep_collinear(mut self, keep: bool) -> Self {
        self.keep_collinear = keep;
        self
    }

    /// Always processes all rings together, never splitting them into independent clusters
    /// (builder form). The output is the same either way; for tests and benchmarks.
    #[doc(hidden)]
    pub fn monolithic(mut self, yes: bool) -> Self {
        self.clustering = if yes {
            Clustering::Never
        } else {
            Clustering::Auto
        };
        self
    }

    /// Splits the rings into independent clusters whenever possible, even when that does not
    /// look worthwhile (builder form). The output is the same either way; for tests.
    #[doc(hidden)]
    pub fn force_clusters(mut self, yes: bool) -> Self {
        self.clustering = if yes {
            Clustering::Always
        } else {
            Clustering::Auto
        };
        self
    }

    /// Adds subject rings and sets the subject fill rule.
    pub fn add_subject(&mut self, rings: &(impl RingSource + ?Sized), rule: FillRule) -> &mut Self {
        self.fill[0] = rule;
        self.add_rings(rings, 0);
        self
    }

    /// Adds clip rings and sets the clip fill rule.
    pub fn add_clip(&mut self, rings: &(impl RingSource + ?Sized), rule: FillRule) -> &mut Self {
        self.fill[1] = rule;
        self.add_rings(rings, 1);
        self
    }

    /// Sets the operation.
    pub fn set_op(&mut self, op: Op) -> &mut Self {
        self.op = op;
        self
    }

    /// Sets whether collinear vertices are kept.
    pub fn set_keep_collinear(&mut self, keep: bool) -> &mut Self {
        self.keep_collinear = keep;
        self
    }

    fn add_rings(&mut self, rings: &(impl RingSource + ?Sized), operand: u8) {
        let edges = &mut self.edges;
        let starts = &mut self.ring_starts;
        let error = &mut self.error;
        rings.visit_rings(&mut |pts, tags| {
            let n = pts.len();
            starts.push(edges.len() as u32);
            if error.is_none()
                && let Some(&p) = pts.iter().find(|p| !p.in_range())
            {
                *error = Some(Error::CoordinateOutOfRange(p));
            }
            for i in 0..n {
                let a = pts[i];
                let b = pts[(i + 1) % n];
                if a != b {
                    let tag = tags.and_then(|t| t.get(i)).copied().unwrap_or(0);
                    edges.push(InEdge { a, b, tag, operand });
                }
            }
        });
    }

    /// Runs the operation, returning the canonical polygon set.
    pub fn execute(&mut self) -> Result<PolygonSet> {
        Ok(self.execute_tree()?.to_polygon_set())
    }

    /// Runs the operation, returning polygons with edge tags.
    pub fn execute_tagged(&mut self) -> Result<Vec<TaggedPolygon>> {
        Ok(self.execute_tree()?.to_tagged_polygons())
    }

    /// Runs the operation, returning the full nesting tree.
    pub fn execute_tree(&mut self) -> Result<PolyTree> {
        if let Some(e) = &self.error {
            return Err(e.clone());
        }
        if self.edges.len() > u32::MAX as usize / 2 {
            return Err(Error::TooLarge);
        }
        let (fill, op) = (self.fill, self.op);
        let inside = move |w: [i32; 2]| op.apply(fill[0].is_inside(w[0]), fill[1].is_inside(w[1]));
        let engine = crate::arrangement::engine();
        if engine != Engine::Auto {
            let boundary = crate::arrangement::boundary(&self.edges, inside, engine);
            return Ok(assemble(boundary, self.keep_collinear));
        }
        if self.clustering != Clustering::Never
            && !ALWAYS_MONOLITHIC.with(|m| m.get())
            && let Some(tree) = crate::cluster::execute(
                &self.edges,
                &self.ring_starts,
                inside,
                self.keep_collinear,
                self.clustering == Clustering::Always,
            )
        {
            return Ok(tree);
        }
        let boundary = compute(&self.edges, inside);
        Ok(assemble(boundary, self.keep_collinear))
    }
}

std::thread_local! {
    static ALWAYS_MONOLITHIC: core::cell::Cell<bool> = const { core::cell::Cell::new(false) };
}

/// Makes every boolean started on the current thread (including those inside offsets and
/// other operations), and every [`clip_paths`], run in one piece, as with
/// [`Boolean::monolithic`]. The output is the
/// same either way; for tests and benchmarks.
#[doc(hidden)]
pub fn set_always_monolithic(yes: bool) {
    ALWAYS_MONOLITHIC.with(|m| m.set(yes));
}

/// Whether rings are split into independent clusters (see `cluster.rs`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
enum Clustering {
    /// When it looks worthwhile.
    #[default]
    Auto,
    Never,
    /// Whenever possible (tests).
    Always,
}

/// Computes the directed boundary edges of the region `{ inside(winding) }`, in sweep
/// order, each with the connected regions on its two sides (see
/// [`boundary`](crate::arrangement::boundary)).
pub(crate) fn compute(edges: &[InEdge], inside: impl Fn([i32; 2]) -> bool + Sync) -> Vec<DirEdge> {
    crate::arrangement::boundary(edges, inside, Engine::Auto)
}

/// Computes `subject op clip` with one fill rule for both operands.
///
/// ```
/// use polyclip::{boolean, FillRule, Op, Ring};
/// let a = Ring::from([(0, 0), (10, 0), (10, 10), (0, 10)]);
/// let b = Ring::from([(5, 5), (15, 5), (15, 15), (5, 15)]);
/// let i = boolean(Op::Intersection, &a, &b, FillRule::NonZero).unwrap();
/// assert_eq!(i[0].outer.signed_area2(), 2 * 25);
/// ```
pub fn boolean(
    op: Op,
    subject: &(impl RingSource + ?Sized),
    clip: &(impl RingSource + ?Sized),
    rule: FillRule,
) -> Result<PolygonSet> {
    Boolean::new()
        .subject(subject, rule)
        .clip(clip, rule)
        .op(op)
        .execute()
}

/// N-ary union of all rings in one pass, under `rule`. Also normalizes arbitrary
/// (self-intersecting, overlapping, any orientation) rings into a canonical polygon set.
pub fn union_all(rings: &(impl RingSource + ?Sized), rule: FillRule) -> Result<PolygonSet> {
    Boolean::new().subject(rings, rule).execute()
}

/// Result of clipping open paths by polygons.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ClippedPaths {
    /// Pieces inside the clip region (including pieces running along its boundary).
    pub inside: Vec<TaggedPath>,
    /// Pieces outside the clip region.
    pub outside: Vec<TaggedPath>,
}

impl ClippedPaths {
    /// Inside pieces without tags.
    pub fn inside_paths(&self) -> Vec<Path> {
        self.inside.iter().map(|p| Path(p.points.clone())).collect()
    }

    /// Outside pieces without tags.
    pub fn outside_paths(&self) -> Vec<Path> {
        self.outside
            .iter()
            .map(|p| Path(p.points.clone()))
            .collect()
    }
}

/// Anything that can be clipped as a set of open paths.
pub trait PathSource {
    /// Calls `f(points, tags)` for every path. `tags`, when present, has one tag per edge.
    fn visit_paths(&self, f: &mut VertexVisitor<'_>);
}

impl PathSource for Path {
    fn visit_paths(&self, f: &mut VertexVisitor<'_>) {
        f(&self.0, None)
    }
}

impl PathSource for TaggedPath {
    fn visit_paths(&self, f: &mut VertexVisitor<'_>) {
        f(&self.points, Some(&self.tags))
    }
}

impl<T: PathSource> PathSource for [T] {
    fn visit_paths(&self, f: &mut VertexVisitor<'_>) {
        for x in self {
            x.visit_paths(f)
        }
    }
}

impl<T: PathSource> PathSource for Vec<T> {
    fn visit_paths(&self, f: &mut VertexVisitor<'_>) {
        self.as_slice().visit_paths(f)
    }
}

impl<T: PathSource + ?Sized> PathSource for &T {
    fn visit_paths(&self, f: &mut VertexVisitor<'_>) {
        (**self).visit_paths(f)
    }
}

/// Clips open paths against the region of `clip` (under `rule`), splitting them into the
/// parts inside and outside.
///
/// Paths are snap-rounded together with the clip edges (so their vertices may move by up to
/// `sqrt(2)/2` units at intersections). Parts running exactly along the clip boundary count
/// as inside. Output pieces keep the direction of their source path and come in source
/// order; consecutive vertices are never equal. Collinear vertices are kept. Each edge keeps
/// the tag of its source edge.
///
/// Only the clip rings near the paths (and the rings near those, transitively) are noded;
/// the others are located once, so a few tracks clipped against a large zone fill cost
/// little more than reading it. The result is identical to noding everything at once.
///
/// ```
/// use polyclip::{clip_paths, FillRule, Path, Ring};
/// let square = Ring::from([(0, 0), (10, 0), (10, 10), (0, 10)]);
/// let line = Path::from([(-5, 5), (15, 5)]);
/// let r = clip_paths(&line, &square, FillRule::NonZero).unwrap();
/// assert_eq!(r.inside_paths(), vec![Path::from([(0, 5), (10, 5)])]);
/// assert_eq!(r.outside.len(), 2);
/// ```
pub fn clip_paths(
    paths: &(impl PathSource + ?Sized),
    clip: &(impl RingSource + ?Sized),
    rule: FillRule,
) -> Result<ClippedPaths> {
    // Clip edges (as `Boolean::add_clip` lists them) and the first edge of every ring.
    let mut n = 0usize;
    clip.visit_rings(&mut |pts, _| n += pts.len());
    let mut edges: Vec<InEdge> = Vec::with_capacity(n);
    let mut ring_starts: Vec<u32> = Vec::new();
    let mut err = None;
    clip.visit_rings(&mut |pts, tags| {
        ring_starts.push(edges.len() as u32);
        if err.is_none()
            && let Some(&p) = pts.iter().find(|p| !p.in_range())
        {
            err = Some(Error::CoordinateOutOfRange(p));
        }
        let tag = |i: usize| tags.and_then(|t| t.get(i)).copied().unwrap_or(0);
        let mut push = |i: usize, a: Point, b: Point| {
            if a != b {
                edges.push(InEdge {
                    a,
                    b,
                    tag: tag(i),
                    operand: 1,
                });
            }
        };
        for (i, w) in pts.windows(2).enumerate() {
            push(i, w[0], w[1]);
        }
        if let (Some(&a), Some(&b)) = (pts.last(), pts.first()) {
            push(pts.len() - 1, a, b);
        }
    });
    if let Some(e) = err {
        return Err(e);
    }
    // Open path edges, remembering where each path starts.
    let mut path_starts: Vec<usize> = Vec::new();
    let mut err = None;
    paths.visit_paths(&mut |pts, tags| {
        path_starts.push(edges.len());
        if err.is_none()
            && let Some(&p) = pts.iter().find(|p| !p.in_range())
        {
            err = Some(Error::CoordinateOutOfRange(p));
        }
        for (i, w) in pts.windows(2).enumerate() {
            if w[0] != w[1] {
                let tag = tags.and_then(|t| t.get(i)).copied().unwrap_or(0);
                edges.push(InEdge {
                    a: w[0],
                    b: w[1],
                    tag,
                    operand: 2,
                });
            }
        }
    });
    if let Some(e) = err {
        return Err(e);
    }
    let first_open = path_starts.first().copied().unwrap_or(edges.len());
    let frags = (!ALWAYS_MONOLITHIC.with(|m| m.get()))
        .then(|| crate::clip::clustered(&edges, &ring_starts, first_open, &path_starts, rule))
        .flatten()
        .unwrap_or_else(|| clip_in_one_piece(&edges, rule));
    // Reassemble pieces per source path.
    let mut res = ClippedPaths::default();
    let mut cur: Option<(bool, TaggedPath)> = None;
    let mut path_idx = 0usize;
    let mut cur_path = usize::MAX;
    for &(a, b, src, ins) in &frags {
        let o = (a, b, edges[src as usize].tag, src);
        while path_idx < path_starts.len() && path_starts[path_idx] <= o.3 as usize {
            path_idx += 1;
        }
        let pid = path_idx - 1;
        let cont = matches!(&cur, Some((f, p)) if *f == ins && cur_path == pid && p.points.last() == Some(&o.0));
        if !cont {
            if let Some((f, p)) = cur.take() {
                if f {
                    res.inside.push(p)
                } else {
                    res.outside.push(p)
                }
            }
            cur = Some((
                ins,
                TaggedPath {
                    points: vec![o.0],
                    tags: Vec::new(),
                },
            ));
            cur_path = pid;
        }
        let p = &mut cur.as_mut().unwrap().1;
        p.points.push(o.1);
        p.tags.push(o.2);
    }
    if let Some((f, p)) = cur.take() {
        if f {
            res.inside.push(p)
        } else {
            res.outside.push(p)
        }
    }
    Ok(res)
}

/// The fragments of the open edges of `edges` with their inside flags (see
/// [`crate::clip::clustered`]), from the arrangement of all edges.
fn clip_in_one_piece(edges: &[InEdge], rule: FillRule) -> Vec<crate::clip::OpenFrag> {
    let Ok(arr) = Arrangement::build(edges, Noding::Snap) else {
        unreachable!("snap rounding never fails")
    };
    let mut inside_flag = vec![false; arr.open_frags.len()];
    for (k, e) in arr.edges.iter().enumerate() {
        if let Some(o) = e.open {
            let (wb, wa) = arr.sides(k);
            inside_flag[o as usize] = rule.is_inside(wb[1]) || rule.is_inside(wa[1]);
        }
    }
    arr.open_frags
        .iter()
        .zip(inside_flag)
        .map(|(f, ins)| (f.a, f.b, f.src, ins))
        .collect()
}
