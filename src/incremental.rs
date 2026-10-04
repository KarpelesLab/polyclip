//! Incremental zone refill: `zone − ⋃ obstacles`, recomputed locally after each change.
//!
//! See [`ZoneFill`] for the API and the exactness argument.

use crate::assemble::{DirEdge, RawRing, link_rings, remove_collinear, rotate_to_min};
use crate::boolean::{FillRule, RingSource};
use crate::error::{Error, Result};
use crate::geom::{
    Point, PolyNode, PolyTree, Polygon, PolygonSet, Rect, Ring, TaggedPolygon, TaggedRing,
};
use crate::node::{Rel, relation, rounded_crossing};
use crate::predicates::{
    cmp_angle, cmp_dir_halfplane, dist2, dot, floor_div, in_segment_interior, orient,
    segment_meets_rect, segments_intersect, sub,
};
use crate::query::ring_area2;
use crate::sweep::{cmp_sweep_edges, sweep};
use core::cmp::Ordering;
use core::hash::{BuildHasherDefault, Hasher};
use std::collections::{BTreeMap, HashMap};

// ---------------------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------------------

/// Incremental zone-fill engine: maintains `zone − ⋃ obstacles` while obstacles are
/// inserted, removed and moved, recomputing only the regions a change can affect.
///
/// The result after any sequence of edits is **identical** (`==`) to recomputing from
/// scratch with
///
/// ```text
/// boolean(Op::Difference, zone, all_obstacles, rule)
/// ```
///
/// (equivalently [`Boolean`](crate::Boolean) with the same fill rule on both operands): the
/// same snap-rounded arrangement, the same canonical rings, tags and nesting. Obstacles are
/// identified by caller-chosen `u64` ids; each obstacle is any [`RingSource`] (rings,
/// polygons, tagged rings, trees, ...), with edge tags carried through as in
/// [`Boolean::execute_tagged`](crate::Boolean::execute_tagged).
///
/// ```
/// use polyclip::*;
///
/// let zone = Ring::from([(0, 0), (1000, 0), (1000, 1000), (0, 1000)]);
/// let sq = |x: i64, y: i64| Ring::from([(x, y), (x + 50, y), (x + 50, y + 50), (x, y + 50)]);
///
/// let mut fill = ZoneFill::new(&zone, FillRule::NonZero).unwrap();
/// fill.insert(1, &sq(100, 100)).unwrap();
/// fill.insert(2, &sq(300, 700)).unwrap();
/// assert_eq!(fill.fill()[0].holes.len(), 2);
///
/// // Move obstacle 2 and remove obstacle 1: only the touched regions are recomputed.
/// fill.update(2, &sq(310, 690)).unwrap();
/// fill.remove(1);
/// let expected = boolean(Op::Difference, &zone, &sq(310, 690), FillRule::NonZero).unwrap();
/// assert_eq!(fill.fill(), expected);
/// ```
///
/// # Editing model
///
/// [`insert`](Self::insert), [`update`](Self::update) and [`remove`](Self::remove) validate
/// their input immediately (coordinates within `±`[`MAX_COORD`](crate::MAX_COORD)) and queue
/// the change; nothing is recomputed until the result is requested
/// ([`fill`](Self::fill), [`fill_tree`](Self::fill_tree), [`fill_tagged`](Self::fill_tagged),
/// [`result`](Self::result)) or [`commit`](Self::commit) is called. All changes queued in
/// between are processed together as one batch. Inserting and then removing the same id
/// before a commit costs nothing.
///
/// # How it works
///
/// The engine keeps the whole state of the boolean pipeline and updates each stage locally:
///
/// 1. **Snap rounding.** The from-scratch boolean snap-rounds all input edges with
///    *selective* hot-pixel activation: crossing pixels are hot, a segment meeting a hot
///    pixel (other than at its own endpoints) is *affected* (rerouted), and every candidate
///    pixel (segment endpoint or rounded crossing) an affected segment meets becomes hot in
///    turn, up to the least fixpoint. A segment's fragments are a function of the segment,
///    its *affected* flag and the candidate pixels within distance 1 of it. The engine keeps
///    the candidate pixels (with endpoint and crossing counts), the hot flags and the
///    affected flags. On an edit it computes the connected component(s) of the realized
///    "segment meets hot pixel" graph that touch the removed segments and pixels (this is
///    where a rounding cascade could have come from), clears them, and re-runs the
///    propagation from the remaining crossing pixels inside them and from everything the
///    new segments introduce. Components untouched by the edit are closed sets generated
///    by crossings that still exist, so they are part of the new least fixpoint unchanged;
///    the recomputed part is the least fixpoint of the rest. The result is exactly the
///    fixpoint of a from-scratch run, however far a cascade propagates. Fragments are then
///    recomputed for every segment whose flag changed or that lies within distance 1 of a
///    pixel that appeared or disappeared, using the same predicates as the batch noder.
/// 2. **Arrangement.** Fragments are merged by exact endpoints into arrangement edges
///    (winding deltas per operand and the provenance tag, with the batch rule). The set `D`
///    of *dirty rectangles* is the bounding box of every inserted or removed ring plus the
///    bounding box of the fragments that changed for each surviving segment. The winding
///    numbers can only change inside `D`: the old and new fragment sets differ by a sum of
///    closed cycles (whole rings, or the old and new fragment chains of one segment, which
///    share their endpoints), each lying inside one rectangle of `D`, and a cycle has
///    winding number zero outside its bounding box. So every arrangement edge not meeting
///    `D` keeps its existence, tag and windings. The edges meeting `D` get new winding
///    numbers in sweep order, each by shooting a ray down to the arrangement edge just
///    below it (exactly the edge the batch sweep would see below it), whose winding is
///    either already recomputed or unchanged.
/// 3. **Rings.** Output rings are cycles of boundary edges, linked at each vertex by a
///    rule that depends only on the boundary edges at that vertex. Rings with an edge
///    meeting `D` or a vertex where a boundary edge changed are dissolved; their unchanged
///    edges plus the new boundary edges meeting `D` are relinked with the batch linker,
///    and canonicalized (collinear vertices, start vertex) as in the batch assembler.
///    Every other ring is unchanged.
/// 4. **Nesting.** Each ring remembers the boundary edge hit by a downward ray from its
///    lowest-leftmost vertex (the batch nesting rule; rays are kept in a spatial index).
///    Only new rings and rings whose ray meets a changed boundary edge are re-shot. A ring's
///    parent follows from the ring owning that edge, sometimes inheriting that ring's
///    parent; parents are re-derived bottom-up for the re-shot rings, the rings whose hit
///    edge now belongs to a new ring, and (transitively) the rings inheriting a parent
///    that changed.
///
/// Every step reproduces the batch definition exactly rather than approximating it, so no
/// tolerance or "margin" is involved, and no verification against a full recompute is
/// needed. (Two batch details matter for exactness and are reproduced: rings cut from the
/// same closed boundary walk at pinch vertices are dissolved together, and relinked edges
/// are walked in sweep order, because splitting a walk at repeated vertices depends on
/// where it starts.) As a safety net, cheap internal consistency checks (balanced boundary
/// degrees at every relinked vertex, a strictly bottom-up nesting order) run on every
/// update, and a failure triggers a full rebuild; [`verify`](Self::verify) compares the
/// whole internal state with a from-scratch recomputation (used by the tests). Very large batches (replacing more than about three quarters of all edges)
/// and large changes in the number of edges also use a full rebuild, which is cheaper there
/// (see [`set_auto_rebuild`](Self::set_auto_rebuild)).
///
/// # Cost
///
/// An edit costs time roughly proportional to the size of the edited obstacles, their
/// neighbourhood, the rounding clusters they touch and the output rings passing through the
/// dirty rectangles (a ring is relinked as a whole). On the 100 mm zone with 5 000
/// obstacles, moving one obstacle takes about 0.1 ms against about 50 ms for a full
/// recompute. Materializing the result ([`fill`](Self::fill) and friends) is linear in the
/// output size; [`result`](Self::result) borrows a cached tree instead, rebuilt after an
/// edit by moving the unchanged rings over from the previous one (no vertex is copied
/// twice).
///
/// The engine is `Send + Sync`, deterministic, and never panics.
#[derive(Clone, Debug)]
pub struct ZoneFill {
    rule: FillRule,
    zone: Vec<InRing>,
    obstacles: BTreeMap<u64, Obstacle>,
    pending: BTreeMap<u64, Option<Vec<InRing>>>,
    st: State,
    /// Total number of input segments (zone and obstacles), committed.
    nseg: usize,
    /// Number of segments when the spatial index was sized.
    nseg_built: usize,
    /// Segments in pending insertions (an upper bound on the next commit's growth).
    pending_segs: usize,
    tree: Option<PolyTree>,
    /// The ring slot and stamp of every node of `tree`.
    tree_stamps: Vec<(u32, u64)>,
    /// The last tree with its ring slots and stamps, after a change: its unchanged rings
    /// are moved into the next tree instead of being copied again.
    stale: Option<(PolyTree, Vec<(u32, u64)>)>,
    rebuilds: u64,
    auto_rebuild: bool,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct InRing {
    pts: Vec<Point>,
    tags: Vec<u64>,
}

#[derive(Clone, Debug)]
struct Obstacle {
    rings: Vec<InRing>,
    segs: Vec<u32>,
}

fn collect_rings(src: &(impl RingSource + ?Sized)) -> Result<(Vec<InRing>, usize)> {
    let mut out = Vec::new();
    let mut err = None;
    let mut n = 0usize;
    src.visit_rings(&mut |pts, tags| {
        if err.is_none()
            && let Some(&p) = pts.iter().find(|p| !p.in_range())
        {
            err = Some(Error::CoordinateOutOfRange(p));
        }
        let t: Vec<u64> = (0..pts.len())
            .map(|i| tags.and_then(|t| t.get(i)).copied().unwrap_or(0))
            .collect();
        n += pts.len();
        out.push(InRing {
            pts: pts.to_vec(),
            tags: t,
        });
    });
    match err {
        Some(e) => Err(e),
        None => Ok((out, n)),
    }
}

/// Upper bound on the number of input segments the engine accepts.
const MAX_SEGS: usize = (u32::MAX / 4) as usize;

impl ZoneFill {
    /// Creates an engine for `zone` with no obstacles. `rule` is the fill rule of both the
    /// zone and the obstacles (as in [`boolean`](crate::boolean)).
    ///
    /// Fails with [`Error::CoordinateOutOfRange`] when a zone coordinate is out of range.
    pub fn new(zone: &(impl RingSource + ?Sized), rule: FillRule) -> Result<Self> {
        let (zone, n) = collect_rings(zone)?;
        if n > MAX_SEGS {
            return Err(Error::TooLarge);
        }
        let mut z = ZoneFill {
            rule,
            zone,
            obstacles: BTreeMap::new(),
            pending: BTreeMap::new(),
            st: State::empty(rule, 0, None),
            nseg: 0,
            nseg_built: 0,
            pending_segs: 0,
            tree: None,
            tree_stamps: Vec::new(),
            stale: None,
            rebuilds: 0,
            auto_rebuild: true,
        };
        z.rebuild();
        Ok(z)
    }

    /// The fill rule of the zone and the obstacles.
    pub fn fill_rule(&self) -> FillRule {
        self.rule
    }

    /// Replaces the zone outline. This rebuilds the whole state (the zone usually touches
    /// everything); obstacles are kept.
    pub fn set_zone(&mut self, zone: &(impl RingSource + ?Sized)) -> Result<()> {
        let (zone, n) = collect_rings(zone)?;
        if n > MAX_SEGS {
            return Err(Error::TooLarge);
        }
        self.apply_pending_to_inputs();
        self.zone = zone;
        self.rebuild();
        Ok(())
    }

    /// Inserts obstacle `id`, replacing any obstacle with the same id. Returns `true` when
    /// an obstacle with this id existed (counting queued changes).
    ///
    /// Fails (leaving the engine unchanged) with [`Error::CoordinateOutOfRange`] when a
    /// coordinate is out of range, or [`Error::TooLarge`] when the total number of edges
    /// would exceed the engine's capacity (about a billion).
    pub fn insert(&mut self, id: u64, obstacle: &(impl RingSource + ?Sized)) -> Result<bool> {
        let (rings, n) = collect_rings(obstacle)?;
        if self.nseg + self.pending_segs + n > MAX_SEGS {
            return Err(Error::TooLarge);
        }
        let existed = self.contains(id);
        self.pending_segs += n;
        self.pending.insert(id, Some(rings));
        Ok(existed)
    }

    /// Replaces obstacle `id` if it exists (counting queued changes); otherwise does
    /// nothing. Returns whether the obstacle existed.
    ///
    /// Fails (leaving the engine unchanged) like [`insert`](Self::insert).
    pub fn update(&mut self, id: u64, obstacle: &(impl RingSource + ?Sized)) -> Result<bool> {
        if !self.contains(id) {
            // Still validate, so that update and insert report the same errors.
            collect_rings(obstacle)?;
            return Ok(false);
        }
        self.insert(id, obstacle)
    }

    /// Removes obstacle `id`. Returns whether it existed (counting queued changes).
    pub fn remove(&mut self, id: u64) -> bool {
        let existed = self.contains(id);
        if existed {
            if self.obstacles.contains_key(&id) {
                self.pending.insert(id, None);
            } else {
                self.pending.remove(&id);
            }
        }
        existed
    }

    /// Removes all obstacles.
    pub fn clear(&mut self) {
        self.pending.clear();
        self.pending_segs = 0;
        for &id in self.obstacles.keys() {
            self.pending.insert(id, None);
        }
    }

    /// Whether obstacle `id` exists (counting queued changes).
    pub fn contains(&self, id: u64) -> bool {
        match self.pending.get(&id) {
            Some(p) => p.is_some(),
            None => self.obstacles.contains_key(&id),
        }
    }

    /// Number of obstacles (counting queued changes).
    pub fn len(&self) -> usize {
        let mut n = self.obstacles.len();
        for (id, p) in &self.pending {
            match (self.obstacles.contains_key(id), p.is_some()) {
                (true, false) => n -= 1,
                (false, true) => n += 1,
                _ => {}
            }
        }
        n
    }

    /// `true` when there are no obstacles (counting queued changes).
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Ids of all obstacles (counting queued changes), in increasing order.
    pub fn ids(&self) -> Vec<u64> {
        let mut v: Vec<u64> = self
            .obstacles
            .keys()
            .filter(|id| !matches!(self.pending.get(id), Some(None)))
            .copied()
            .collect();
        v.extend(
            self.pending
                .iter()
                .filter(|(id, p)| p.is_some() && !self.obstacles.contains_key(id))
                .map(|(id, _)| *id),
        );
        v.sort_unstable();
        v
    }

    /// Enables or disables automatic full rebuilds (enabled by default).
    ///
    /// When enabled, a batch replacing more than about three quarters of all edges, or a commit
    /// that doubles (or quarters) the total number of edges since the spatial index was last
    /// sized, is processed by rebuilding the whole state, which is
    /// faster in those cases. When disabled, every commit takes the incremental path (the
    /// result is identical either way; only speed differs). Rebuilds after a failed
    /// internal consistency check happen regardless.
    pub fn set_auto_rebuild(&mut self, enabled: bool) {
        self.auto_rebuild = enabled;
    }

    /// Number of full rebuilds performed so far (initial construction, zone changes, large
    /// batches and consistency fallbacks). Diagnostic only.
    pub fn rebuild_count(&self) -> u64 {
        self.rebuilds
    }

    /// Commits queued changes, then checks the engine's internal state against a
    /// from-scratch recomputation of every stage (snap-rounding fixpoint and fragments,
    /// arrangement windings, rings and nesting). Returns a description of the first
    /// difference found.
    ///
    /// This is slow (a full rebuild and more); it exists for tests and debugging. The
    /// result returned by [`fill`](Self::fill) is always exact regardless.
    pub fn verify(&mut self) -> core::result::Result<(), String> {
        self.commit();
        if let Some(m) = self.st.check_snap() {
            return Err(format!("snap rounding: {m}"));
        }
        if let Some(m) = self.st.check_windings() {
            return Err(format!("windings: {m}"));
        }
        if let Some(m) = self.st.check_rings() {
            return Err(format!("rings: {m}"));
        }
        Ok(())
    }

    /// Processes all queued changes.
    pub fn commit(&mut self) {
        if self.pending.is_empty() {
            return;
        }
        // Drop no-op changes (re-inserting identical geometry).
        let pending = core::mem::take(&mut self.pending);
        self.pending_segs = 0;
        let mut changes: Vec<(u64, Option<Vec<InRing>>)> = Vec::new();
        let mut add_n = 0usize;
        let mut rem_n = 0usize;
        for (id, p) in pending {
            let old = self.obstacles.get(&id);
            if let (Some(o), Some(n)) = (old, &p)
                && o.rings == *n
            {
                continue;
            }
            if old.is_none() && p.is_none() {
                continue;
            }
            rem_n += old.map_or(0, |o| o.segs.len());
            add_n += p
                .as_ref()
                .map_or(0, |r| r.iter().map(|r| r.pts.len()).sum());
            changes.push((id, p));
        }
        if changes.is_empty() {
            return;
        }
        self.invalidate_tree();
        let new_total = self.nseg + add_n - rem_n.min(self.nseg + add_n);
        let big = self.auto_rebuild
            && ((add_n + rem_n) * 4 > self.nseg.max(1) * 3
                || new_total > 2 * self.nseg_built + 1024
                || new_total * 4 + 1024 < self.nseg_built);
        if big {
            for (id, p) in changes {
                match p {
                    Some(rings) => {
                        self.obstacles.insert(
                            id,
                            Obstacle {
                                rings,
                                segs: Vec::new(),
                            },
                        );
                    }
                    None => {
                        self.obstacles.remove(&id);
                    }
                }
            }
            self.rebuild();
            return;
        }
        // Incremental path.
        let mut rem_segs: Vec<u32> = Vec::new();
        let mut rects: Vec<Rect> = Vec::new();
        let mut add: Vec<NewSeg> = Vec::new();
        let mut add_owner: Vec<(u64, usize)> = Vec::new();
        for (id, p) in &changes {
            if let Some(o) = self.obstacles.get(id) {
                rem_segs.extend_from_slice(&o.segs);
                for r in &o.rings {
                    if let Some(bb) = Rect::of_points(r.pts.iter()) {
                        rects.push(bb);
                    }
                }
            }
            if let Some(rings) = p {
                let start = add.len();
                for r in rings {
                    push_ring_segs(r, 1, &mut add);
                    if let Some(bb) = Rect::of_points(r.pts.iter()) {
                        rects.push(bb);
                    }
                }
                add_owner.push((*id, add.len() - start));
            }
        }
        let ok = self.st.apply(&rem_segs, &add, rects);
        // Commit the inputs.
        let mut new_ids = match &ok {
            Ok(ids) => ids.clone(),
            Err(_) => Vec::new(),
        }
        .into_iter();
        let mut owner = add_owner.into_iter();
        for (id, p) in changes {
            match p {
                Some(rings) => {
                    let n = owner.next().map_or(0, |o| o.1);
                    let segs: Vec<u32> = new_ids.by_ref().take(n).collect();
                    self.obstacles.insert(id, Obstacle { rings, segs });
                }
                None => {
                    self.obstacles.remove(&id);
                }
            }
        }
        self.nseg = self.st.live_segs;
        if ok.is_err() {
            self.rebuild();
        }
    }

    /// Commits queued changes and returns the current fill as a canonical polygon set,
    /// identical to `boolean(Op::Difference, zone, obstacles, rule)`.
    pub fn fill(&mut self) -> PolygonSet {
        self.commit();
        if let Some(t) = &self.tree {
            return t.to_polygon_set();
        }
        self.st
            .flatten()
            .into_iter()
            .map(|(o, hs)| Polygon {
                outer: Ring(o.pts.clone()),
                holes: hs.iter().map(|h| Ring(h.pts.clone())).collect(),
            })
            .collect()
    }

    /// Commits queued changes and returns the current fill as a nesting tree, identical to
    /// [`Boolean::execute_tree`](crate::Boolean::execute_tree) on the same input.
    pub fn fill_tree(&mut self) -> PolyTree {
        self.result().clone()
    }

    /// Commits queued changes and returns the current fill with edge tags, identical to
    /// [`Boolean::execute_tagged`](crate::Boolean::execute_tagged) on the same input.
    pub fn fill_tagged(&mut self) -> Vec<TaggedPolygon> {
        self.commit();
        if let Some(t) = &self.tree {
            return t.to_tagged_polygons();
        }
        let tr = |r: &RawRing| TaggedRing {
            points: r.pts.clone(),
            tags: r.tags.clone(),
        };
        self.st
            .flatten()
            .into_iter()
            .map(|(o, hs)| TaggedPolygon {
                outer: tr(o),
                holes: hs.into_iter().map(tr).collect(),
            })
            .collect()
    }

    /// Commits queued changes and borrows the current fill tree (cached until the next
    /// change).
    pub fn result(&mut self) -> &PolyTree {
        self.commit();
        if self.tree.is_none() {
            let (t, stamps) = self.st.tree(self.stale.take());
            self.tree = Some(t);
            self.tree_stamps = stamps;
        }
        self.tree.get_or_insert_with(PolyTree::default)
    }

    /// Drops the cached tree, keeping it (and its ring stamps) for reuse.
    fn invalidate_tree(&mut self) {
        if let Some(t) = self.tree.take() {
            self.stale = Some((t, core::mem::take(&mut self.tree_stamps)));
        }
    }

    /// Applies queued changes to the stored inputs without computing anything.
    fn apply_pending_to_inputs(&mut self) {
        for (id, p) in core::mem::take(&mut self.pending) {
            match p {
                Some(rings) => {
                    self.obstacles.insert(
                        id,
                        Obstacle {
                            rings,
                            segs: Vec::new(),
                        },
                    );
                }
                None => {
                    self.obstacles.remove(&id);
                }
            }
        }
        self.pending_segs = 0;
    }

    /// Rebuilds the whole state from the stored inputs.
    fn rebuild(&mut self) {
        self.apply_pending_to_inputs();
        self.invalidate_tree();
        self.rebuilds += 1;
        let mut segs: Vec<NewSeg> = Vec::new();
        for r in &self.zone {
            push_ring_segs(r, 0, &mut segs);
        }
        let nz = segs.len();
        let mut counts = Vec::with_capacity(self.obstacles.len());
        for o in self.obstacles.values() {
            let s = segs.len();
            for r in &o.rings {
                push_ring_segs(r, 1, &mut segs);
            }
            counts.push(segs.len() - s);
        }
        let ids = State::build(self.rule, &segs, &mut self.st);
        let mut it = ids.into_iter().skip(nz);
        for (o, n) in self.obstacles.values_mut().zip(counts) {
            o.segs = it.by_ref().take(n).collect();
        }
        self.nseg = self.st.live_segs;
        self.nseg_built = self.nseg;
    }
}

/// A segment to insert.
#[derive(Clone, Copy, Debug)]
struct NewSeg {
    a: Point,
    b: Point,
    tag: u64,
    operand: u8,
}

fn push_ring_segs(r: &InRing, operand: u8, out: &mut Vec<NewSeg>) {
    let n = r.pts.len();
    for i in 0..n {
        let a = r.pts[i];
        let b = r.pts[(i + 1) % n];
        if a != b {
            out.push(NewSeg {
                a,
                b,
                tag: r.tags.get(i).copied().unwrap_or(0),
                operand,
            });
        }
    }
}

// ---------------------------------------------------------------------------------------
// Hashing and spatial index
// ---------------------------------------------------------------------------------------

/// Small deterministic multiplicative hasher (no random state: iteration orders are
/// reproducible, although results never depend on them).
#[derive(Clone, Copy, Default)]
struct FxHasher(u64);

impl FxHasher {
    #[inline]
    fn add(&mut self, i: u64) {
        self.0 = (self.0.rotate_left(5) ^ i).wrapping_mul(0x51_7c_c1_b7_27_22_0a_95);
    }
}

impl Hasher for FxHasher {
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.add(b as u64);
        }
    }
    #[inline]
    fn write_u8(&mut self, i: u8) {
        self.add(i as u64)
    }
    #[inline]
    fn write_u32(&mut self, i: u32) {
        self.add(i as u64)
    }
    #[inline]
    fn write_u64(&mut self, i: u64) {
        self.add(i)
    }
    #[inline]
    fn write_i64(&mut self, i: i64) {
        self.add(i as u64)
    }
    #[inline]
    fn write_usize(&mut self, i: usize) {
        self.add(i as u64)
    }
    #[inline]
    fn finish(&self) -> u64 {
        let h = self.0;
        (h ^ (h >> 29)).wrapping_mul(0xbf58_476d_1ce4_e5b9) ^ (h >> 32)
    }
}

type FxMap<K, V> = HashMap<K, V, BuildHasherDefault<FxHasher>>;

type Cell = (i64, i64);

#[inline]
fn cell_of(p: Point, sh: u32) -> Cell {
    (p.x >> sh, p.y >> sh)
}

#[inline]
fn cell_rect_grown(c: Cell, sh: u32) -> Rect {
    Rect {
        min: Point::new((c.0 << sh) - 1, (c.1 << sh) - 1),
        max: Point::new((c.0 + 1) << sh, (c.1 + 1) << sh),
    }
}

#[inline]
fn ceil_div(a: i128, b: i128) -> i128 {
    -floor_div(-a, b)
}

/// Segments listed in more cells than this go to the grid's "long" list instead.
const LONG_CELLS: u64 = 1 << 14;

/// Approximate number of cells [`for_seg_cells`] visits.
fn seg_cell_count(a: Point, b: Point, sh: u32) -> u64 {
    let (x0, x1) = (a.x.min(b.x), a.x.max(b.x));
    let (y0, y1) = (a.y.min(b.y), a.y.max(b.y));
    let nx = (((x1 + 1) >> sh) - ((x0 - 1) >> sh) + 1) as u64;
    let ny = (((y1 + 1) >> sh) - ((y0 - 1) >> sh) + 1) as u64;
    if nx <= 2 || ny <= 2 || a.x == b.x || a.y == b.y {
        nx.saturating_mul(ny)
    } else {
        3 * (nx + ny)
    }
}

/// Calls `f` for every cell whose closed square grown by one unit the segment may meet (a
/// superset of the cells it does meet). Any point within distance 1 of the segment lies in
/// one of these cells.
fn for_seg_cells(a: Point, b: Point, sh: u32, f: &mut impl FnMut(Cell)) {
    let (x0, x1) = (a.x.min(b.x), a.x.max(b.x));
    let (y0, y1) = (a.y.min(b.y), a.y.max(b.y));
    let (cx0, cx1) = ((x0 - 1) >> sh, (x1 + 1) >> sh);
    let (cy0, cy1) = ((y0 - 1) >> sh, (y1 + 1) >> sh);
    if cx1 - cx0 <= 1 || cy1 - cy0 <= 1 || a.x == b.x || a.y == b.y {
        for cy in cy0..=cy1 {
            for cx in cx0..=cx1 {
                f((cx, cy));
            }
        }
        return;
    }
    let (p, q) = if a.x < b.x { (a, b) } else { (b, a) };
    let dx = (q.x - p.x) as i128;
    let dy = (q.y - p.y) as i128;
    for cx in cx0..=cx1 {
        let xa = ((cx << sh) - 1).max(p.x);
        let xb = ((cx + 1) << sh).min(q.x);
        if xa > xb {
            continue;
        }
        let na = (xa - p.x) as i128 * dy;
        let nb = (xb - p.x) as i128 * dy;
        let (nlo, nhi) = if dy >= 0 { (na, nb) } else { (nb, na) };
        let ylo = p.y + floor_div(nlo, dx) as i64 - 1;
        let yhi = p.y + ceil_div(nhi, dx) as i64 + 1;
        for cy in (ylo >> sh).max(cy0)..=(yhi >> sh).min(cy1) {
            f((cx, cy));
        }
    }
}

/// Uniform grid of `u32` ids: a dense array of cells over a fixed extent (chosen when the
/// state is built) plus a hash map for cells outside it. Segments are listed in every cell
/// near them (see [`for_seg_cells`]), points in the cell containing them; segments that
/// would need too many cells live in `long` and are visited by every query.
#[derive(Clone, Debug, Default)]
struct Grid {
    x0: i64,
    y0: i64,
    nx: i64,
    ny: i64,
    dense: Vec<Vec<u32>>,
    sparse: FxMap<Cell, Vec<u32>>,
    long: Vec<u32>,
    /// Number of non-empty cells.
    occupied: usize,
    /// A lower bound of the rows of all cells ever used.
    min_row: i64,
}

/// Largest dense part of a grid, in cells.
const MAX_DENSE: i64 = 1 << 22;

impl Grid {
    /// An empty grid whose dense part covers the cells of `extent` (when not too large).
    fn new(extent: Option<Rect>, sh: u32) -> Self {
        let mut g = Grid {
            min_row: i64::MAX,
            ..Grid::default()
        };
        if let Some(r) = extent {
            let (cx0, cy0) = ((r.min.x >> sh) - 1, (r.min.y >> sh) - 1);
            let (cx1, cy1) = ((r.max.x >> sh) + 1, (r.max.y >> sh) + 1);
            let (nx, ny) = (cx1 - cx0 + 1, cy1 - cy0 + 1);
            if nx.saturating_mul(ny) <= MAX_DENSE {
                g.x0 = cx0;
                g.y0 = cy0;
                g.nx = nx;
                g.ny = ny;
                g.dense = vec![Vec::new(); (nx * ny) as usize];
            }
        }
        g
    }

    #[inline]
    fn slot(&self, c: Cell) -> Option<usize> {
        let (dx, dy) = (c.0 - self.x0, c.1 - self.y0);
        if dx >= 0 && dx < self.nx && dy >= 0 && dy < self.ny {
            Some((dy * self.nx + dx) as usize)
        } else {
            None
        }
    }

    #[inline]
    fn get(&self, c: Cell) -> &[u32] {
        match self.slot(c) {
            Some(i) => &self.dense[i],
            None => self.sparse.get(&c).map_or(&[], |v| v.as_slice()),
        }
    }

    #[inline]
    fn add(&mut self, c: Cell, id: u32) {
        self.min_row = self.min_row.min(c.1);
        let v = match self.slot(c) {
            Some(i) => &mut self.dense[i],
            None => self.sparse.entry(c).or_default(),
        };
        if v.is_empty() {
            self.occupied += 1;
        }
        v.push(id);
    }

    #[inline]
    fn del(&mut self, c: Cell, id: u32) {
        let v = match self.slot(c) {
            Some(i) => Some(&mut self.dense[i]),
            None => self.sparse.get_mut(&c),
        };
        if let Some(v) = v
            && let Some(i) = v.iter().position(|&x| x == id)
        {
            v.swap_remove(i);
            if v.is_empty() {
                self.occupied -= 1;
                if self.slot(c).is_none() {
                    self.sparse.remove(&c);
                }
            }
        }
    }

    /// Calls `f(cell, items)` for every non-empty cell.
    fn for_each_cell(&self, f: &mut impl FnMut(Cell, &[u32])) {
        for (i, v) in self.dense.iter().enumerate() {
            if !v.is_empty() {
                let i = i as i64;
                f((self.x0 + i % self.nx, self.y0 + i / self.nx), v);
            }
        }
        for (&c, v) in &self.sparse {
            f(c, v);
        }
    }

    /// Registers segment `id`; returns `true` when it went to the long list.
    fn add_seg(&mut self, id: u32, a: Point, b: Point, sh: u32) -> bool {
        if seg_cell_count(a, b, sh) > LONG_CELLS {
            self.long.push(id);
            self.min_row = self.min_row.min((a.y.min(b.y) - 1) >> sh);
            return true;
        }
        for_seg_cells(a, b, sh, &mut |c| self.add(c, id));
        false
    }

    fn del_seg(&mut self, id: u32, a: Point, b: Point, sh: u32, long: bool) {
        if long {
            if let Some(i) = self.long.iter().position(|&x| x == id) {
                self.long.swap_remove(i);
            }
            return;
        }
        for_seg_cells(a, b, sh, &mut |c| self.del(c, id));
    }

    /// Calls `f(cell, items)` for every non-empty cell near segment `a-b` (every cell
    /// containing a point within distance 1 of it is visited exactly once). Long-list items
    /// are not included.
    fn near_seg(&self, a: Point, b: Point, sh: u32, f: &mut impl FnMut(Cell, &[u32])) {
        if seg_cell_count(a, b, sh) > LONG_CELLS.min(self.occupied as u64 + 16) {
            self.for_each_cell(&mut |c, v| {
                if segment_meets_rect(a, b, &cell_rect_grown(c, sh)) {
                    f(c, v);
                }
            });
            return;
        }
        for_seg_cells(a, b, sh, &mut |c| {
            let v = self.get(c);
            if !v.is_empty() {
                f(c, v);
            }
        });
    }

    /// Calls `f(items)` for every non-empty cell meeting the closed rectangle `r`.
    fn in_rect(&self, r: &Rect, sh: u32, f: &mut impl FnMut(&[u32])) {
        let (cx0, cx1) = (r.min.x >> sh, r.max.x >> sh);
        let (cy0, cy1) = (r.min.y >> sh, r.max.y >> sh);
        let n = ((cx1 - cx0 + 1) as u128) * ((cy1 - cy0 + 1) as u128);
        if n > self.occupied as u128 + 16 {
            self.for_each_cell(&mut |c, v| {
                if c.0 >= cx0 && c.0 <= cx1 && c.1 >= cy0 && c.1 <= cy1 {
                    f(v);
                }
            });
            return;
        }
        for cy in cy0..=cy1 {
            for cx in cx0..=cx1 {
                let v = self.get((cx, cy));
                if !v.is_empty() {
                    f(v);
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------------------

const NONE: u32 = u32::MAX;

/// An input segment with its snap-rounding state.
#[derive(Clone, Debug)]
struct Seg {
    a: Point,
    b: Point,
    bb: Rect,
    tag: u64,
    operand: u8,
    live: bool,
    long: bool,
    affected: bool,
    /// Fragment chain from `a` to `b` (consecutive vertices distinct).
    chain: Vec<Point>,
}

/// A candidate pixel (segment endpoint and/or rounded proper crossing).
#[derive(Clone, Debug)]
struct Pix {
    p: Point,
    ends: u32,
    cross: u32,
    hot: bool,
    live: bool,
}

/// A fragment contribution to an arrangement edge.
#[derive(Clone, Copy, Debug)]
struct Contrib {
    operand: u8,
    sign: i8,
    tag: u64,
    count: u32,
}

/// An arrangement edge (all coincident fragments merged), `lo < hi`.
#[derive(Clone, Debug)]
struct Edge {
    lo: Point,
    hi: Point,
    contribs: Vec<Contrib>,
    delta: [i32; 2],
    tag: u64,
    /// Winding numbers just below.
    below: [i32; 2],
    /// Boundary direction: 1 when the result lies above (output edge `lo -> hi`), -1 when
    /// below (`hi -> lo`), 0 when not a boundary edge.
    bdir: i8,
    /// Ring of a boundary edge.
    ring: u32,
    live: bool,
    long: bool,
}

impl Edge {
    #[inline]
    fn active(&self) -> bool {
        self.live && self.delta != [0, 0]
    }
    #[inline]
    fn ends(&self) -> (Point, Point) {
        if self.bdir > 0 {
            (self.lo, self.hi)
        } else {
            (self.hi, self.lo)
        }
    }
}

/// What a ring's nesting depends on: the boundary edge just below its query edge.
#[derive(Clone, Copy, Debug)]
struct Wit {
    edge: u32,
    v: Point,
    /// Lower end (rounded down) of the vertical segment from `v` to the hit point.
    ylo: i64,
}

#[derive(Clone, Debug)]
struct RingRec {
    edges: Vec<u32>,
    raw: RawRing,
    is_hole: bool,
    /// `false` for degenerate loops (dropped from the output, as in batch assembly).
    real: bool,
    live: bool,
    wit: Wit,
    query_dir: Point,
    /// Whether the witness ray is registered in `wit_grid` / `wit_of` (and in its long list).
    wit_reg: bool,
    wit_long: bool,
    /// Next ring split from the same closed boundary walk (circular list; itself when the
    /// walk gave one ring). Such rings are always dissolved together.
    gnext: u32,
    /// Unique among the rings ever created by the engine (see `State::next_stamp`).
    stamp: u64,
}

#[derive(Clone, Debug)]
struct State {
    sh: u32,
    rule: FillRule,
    segs: Vec<Seg>,
    seg_free: Vec<u32>,
    seg_grid: Grid,
    live_segs: usize,
    pix: Vec<Pix>,
    pix_free: Vec<u32>,
    pix_map: FxMap<Point, u32>,
    pix_grid: Grid,
    edges: Vec<Edge>,
    edge_free: Vec<u32>,
    edge_map: FxMap<(Point, Point), u32>,
    edge_grid: Grid,
    outdeg: FxMap<Point, u32>,
    rings: Vec<RingRec>,
    ring_free: Vec<u32>,
    /// Next ring stamp: every ring created gets a new one (kept across rebuilds), so a
    /// ring with the same stamp as a node of an earlier tree has the same vertices.
    next_stamp: u64,
    parent: Vec<Option<u32>>,
    /// Witness rays of the rings (vertical segments), by ring id.
    wit_grid: Grid,
    /// Rings whose witness is a given boundary edge.
    wit_of: FxMap<u32, Vec<u32>>,
    // Scratch marks (epoch stamps).
    epoch: u32,
    seg_mark: Vec<u32>,
    pix_mark: Vec<u32>,
    edge_mark: Vec<u32>,
    edge_mark2: Vec<u32>,
    ring_mark: Vec<u32>,
    // Deferred frees (ids are not reused within one update).
    dead_pix: Vec<u32>,
    dead_edges: Vec<u32>,
    dead_segs: Vec<u32>,
}

/// Relation of a segment to a candidate pixel.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum R {
    Far,
    Near,
    Meets,
    Own,
}

#[inline]
fn rel(s: &Seg, p: Point) -> R {
    if p == s.a || p == s.b {
        return R::Own;
    }
    match relation(s.a, s.b, &s.bb, p) {
        Rel::Far => R::Far,
        Rel::Near => R::Near,
        Rel::Meets => R::Meets,
    }
}

/// Cell size (as a power of two) for `n` segments in `bb`: about one cell per segment, and
/// never many more cells than segments (thin, wide boxes).
fn choose_shift(bb: Option<Rect>, n: usize) -> u32 {
    let Some(bb) = bb else { return 16 };
    let w = (bb.max.x - bb.min.x) as f64 + 1.0;
    let h = (bb.max.y - bb.min.y) as f64 + 1.0;
    let s = (w * h / (n.max(1) as f64)).sqrt().max(1.0);
    let s = if s.is_finite() { s as u64 } else { 1 << 40 };
    let mut sh = (63 - s.max(1).leading_zeros()).min(41);
    let cells = |sh: u32| {
        let nx = ((bb.max.x >> sh) - (bb.min.x >> sh) + 3) as u128;
        let ny = ((bb.max.y >> sh) - (bb.min.y >> sh) + 3) as u128;
        nx * ny
    };
    while sh < 41 && cells(sh) > 4 * n as u128 + 1024 {
        sh += 1;
    }
    sh
}

/// Exact `y` of a non-vertical edge at abscissa `x` as a fraction `(num, den)`, `den > 0`.
#[inline]
fn y_at(lo: Point, hi: Point, x: i64) -> (i128, i128) {
    let dx = (hi.x - lo.x) as i128;
    let dy = (hi.y - lo.y) as i128;
    (lo.y as i128 * dx + (x - lo.x) as i128 * dy, dx)
}

impl State {
    fn empty(rule: FillRule, sh: u32, extent: Option<Rect>) -> State {
        State {
            sh,
            rule,
            segs: Vec::new(),
            seg_free: Vec::new(),
            seg_grid: Grid::new(extent, sh),
            live_segs: 0,
            pix: Vec::new(),
            pix_free: Vec::new(),
            pix_map: FxMap::default(),
            pix_grid: Grid::new(extent, sh),
            edges: Vec::new(),
            edge_free: Vec::new(),
            edge_map: FxMap::default(),
            edge_grid: Grid::new(extent, sh),
            outdeg: FxMap::default(),
            rings: Vec::new(),
            ring_free: Vec::new(),
            next_stamp: 0,
            parent: Vec::new(),
            wit_grid: Grid::new(extent, sh),
            wit_of: FxMap::default(),
            epoch: 0,
            seg_mark: Vec::new(),
            pix_mark: Vec::new(),
            edge_mark: Vec::new(),
            edge_mark2: Vec::new(),
            ring_mark: Vec::new(),
            dead_pix: Vec::new(),
            dead_edges: Vec::new(),
            dead_segs: Vec::new(),
        }
    }

    #[inline]
    fn inside(&self, w: [i32; 2]) -> bool {
        self.rule.is_inside(w[0]) && !self.rule.is_inside(w[1])
    }

    fn next_epoch(&mut self) -> u32 {
        self.epoch = self.epoch.wrapping_add(1);
        if self.epoch == 0 {
            for v in [
                &mut self.seg_mark,
                &mut self.pix_mark,
                &mut self.edge_mark,
                &mut self.edge_mark2,
                &mut self.ring_mark,
            ] {
                v.iter_mut().for_each(|x| *x = 0);
            }
            self.epoch = 1;
        }
        self.seg_mark.resize(self.segs.len(), 0);
        self.pix_mark.resize(self.pix.len(), 0);
        self.edge_mark.resize(self.edges.len(), 0);
        self.edge_mark2.resize(self.edges.len(), 0);
        self.ring_mark.resize(self.rings.len(), 0);
        self.epoch
    }

    // ----- batch construction ------------------------------------------------------------

    /// Builds the state for `input` from scratch (into `st`, reusing nothing). Returns the
    /// segment id of every input segment.
    fn build(rule: FillRule, input: &[NewSeg], st: &mut State) -> Vec<u32> {
        let bb = input
            .iter()
            .map(|s| Rect::new(s.a, s.b))
            .reduce(|a, b| a.union(&b));
        let sh = choose_shift(bb, input.len());
        let stamp = st.next_stamp;
        *st = State::empty(rule, sh, bb);
        st.next_stamp = stamp;
        st.segs.reserve(input.len());
        st.pix_map.reserve(input.len() + input.len() / 8);
        st.edge_map.reserve(input.len() + input.len() / 4);
        st.edges.reserve(input.len() + input.len() / 4);
        let mut created = Vec::new();
        let mut new_cross = Vec::new();
        let mut ids = Vec::with_capacity(input.len());
        for s in input {
            let id = st.alloc_seg(s);
            st.register_seg(id, &mut created, &mut new_cross);
            ids.push(id);
        }
        // Fixpoint from all crossing pixels.
        let mut pq: Vec<u32> = Vec::new();
        let mut sq: Vec<u32> = Vec::new();
        for k in 0..st.pix.len() {
            if st.pix[k].live && st.pix[k].cross > 0 && !st.pix[k].hot {
                st.pix[k].hot = true;
                pq.push(k as u32);
            }
        }
        st.propagate(&mut pq, &mut sq, &mut |_, _| {});
        // Fragments and arrangement.
        let mut scratch = Scratch::default();
        for id in 0..st.segs.len() {
            let chain = st.compute_chain(id as u32, &mut scratch);
            let s = &st.segs[id];
            let (op, tag) = (s.operand, s.tag);
            for w in chain.windows(2) {
                st.frag(w[0], w[1], op, tag, true);
            }
            st.segs[id].chain = chain;
        }
        // Windings by one sweep.
        let mut order: Vec<u32> = (0..st.edges.len() as u32)
            .filter(|&e| st.edges[e as usize].active())
            .collect();
        st.sort_sweep(&mut order);
        let segs: Vec<(Point, Point)> = order
            .iter()
            .map(|&e| (st.edges[e as usize].lo, st.edges[e as usize].hi))
            .collect();
        let mut below = vec![[0i32; 2]; order.len()];
        {
            let edges = &st.edges;
            sweep(&segs, |k, b| {
                if let Some(b) = b {
                    let eb = &edges[order[b as usize] as usize];
                    below[k as usize] = [
                        below[b as usize][0] + eb.delta[0],
                        below[b as usize][1] + eb.delta[1],
                    ];
                }
            });
        }
        let mut bnd: Vec<u32> = Vec::new();
        for (k, &e) in order.iter().enumerate() {
            st.set_winding(e, below[k]);
            if st.edges[e as usize].bdir != 0 {
                bnd.push(e);
                let (f, _) = st.edges[e as usize].ends();
                *st.outdeg.entry(f).or_insert(0) += 1;
            }
        }
        let mut new_rings = Vec::new();
        // A fresh build is always balanced.
        let _ = st.link(&bnd, &mut new_rings);
        for &r in &new_rings {
            st.shoot_ring(r);
        }
        if st.compute_parents().is_err() {
            // Cannot happen for a consistent arrangement; keep a valid (flat) nesting.
            st.parent = vec![None; st.rings.len()];
        }
        ids
    }

    // ----- segments and pixels -----------------------------------------------------------

    fn alloc_seg(&mut self, s: &NewSeg) -> u32 {
        let rec = Seg {
            a: s.a,
            b: s.b,
            bb: Rect::new(s.a, s.b),
            tag: s.tag,
            operand: s.operand,
            live: true,
            long: false,
            affected: false,
            chain: Vec::new(),
        };
        self.live_segs += 1;
        if let Some(id) = self.seg_free.pop() {
            self.segs[id as usize] = rec;
            id
        } else {
            self.segs.push(rec);
            self.segs.len() as u32 - 1
        }
    }

    fn pix_id(&mut self, p: Point, created: &mut Vec<Point>) -> u32 {
        if let Some(&k) = self.pix_map.get(&p) {
            return k;
        }
        let rec = Pix {
            p,
            ends: 0,
            cross: 0,
            hot: false,
            live: true,
        };
        let k = if let Some(k) = self.pix_free.pop() {
            self.pix[k as usize] = rec;
            k
        } else {
            self.pix.push(rec);
            self.pix.len() as u32 - 1
        };
        self.pix_map.insert(p, k);
        self.pix_grid.add(cell_of(p, self.sh), k);
        created.push(p);
        k
    }

    fn pix_release(&mut self, k: u32, deleted: &mut Vec<Point>) {
        let px = &mut self.pix[k as usize];
        if px.ends == 0 && px.cross == 0 && px.live {
            px.live = false;
            px.hot = false;
            let p = px.p;
            self.pix_map.remove(&p);
            self.pix_grid.del(cell_of(p, self.sh), k);
            self.dead_pix.push(k);
            deleted.push(p);
        }
    }

    /// Rounded proper crossings of segment `s` with the registered segments, each pair
    /// reported once (in the cell containing the crossing pixel).
    fn crossings(&self, s: u32, out: &mut Vec<Point>) {
        out.clear();
        let sg = &self.segs[s as usize];
        let (a, b, bb) = (sg.a, sg.b, sg.bb);
        let sh = self.sh;
        let test = |t: u32, cell: Option<Cell>, out: &mut Vec<Point>| {
            if t == s {
                return;
            }
            let tg = &self.segs[t as usize];
            if !tg.live || !bb.intersects(&tg.bb) {
                return;
            }
            let (p, q) = (tg.a, tg.b);
            let o1 = orient(a, b, p).signum();
            let o2 = orient(a, b, q).signum();
            if o1 * o2 >= 0 {
                return;
            }
            let o3 = orient(p, q, a).signum();
            let o4 = orient(p, q, b).signum();
            if o3 * o4 < 0 {
                let x = rounded_crossing(a, b, p, q);
                if cell.is_none_or(|c| cell_of(x, sh) == c) {
                    out.push(x);
                }
            }
        };
        self.seg_grid.near_seg(a, b, sh, &mut |c, items| {
            for &t in items {
                test(t, Some(c), out);
            }
        });
        for &t in &self.seg_grid.long {
            test(t, None, out);
        }
    }

    fn register_seg(&mut self, id: u32, created: &mut Vec<Point>, new_cross: &mut Vec<u32>) {
        let (a, b) = (self.segs[id as usize].a, self.segs[id as usize].b);
        let long = self.seg_grid.add_seg(id, a, b, self.sh);
        self.segs[id as usize].long = long;
        for p in [a, b] {
            let k = self.pix_id(p, created);
            self.pix[k as usize].ends += 1;
        }
        let mut xs = Vec::new();
        self.crossings(id, &mut xs);
        for x in xs {
            let k = self.pix_id(x, created);
            let px = &mut self.pix[k as usize];
            if px.cross == 0 {
                new_cross.push(k);
            }
            px.cross += 1;
        }
    }

    fn unregister_seg(&mut self, id: u32, deleted: &mut Vec<Point>) {
        let s = &self.segs[id as usize];
        let (a, b, long) = (s.a, s.b, s.long);
        self.seg_grid.del_seg(id, a, b, self.sh, long);
        self.segs[id as usize].live = false;
        let mut xs = Vec::new();
        self.crossings(id, &mut xs);
        for x in xs {
            if let Some(&k) = self.pix_map.get(&x) {
                let px = &mut self.pix[k as usize];
                px.cross = px.cross.saturating_sub(1);
                self.pix_release(k, deleted);
            }
        }
        for p in [a, b] {
            if let Some(&k) = self.pix_map.get(&p) {
                let px = &mut self.pix[k as usize];
                px.ends = px.ends.saturating_sub(1);
                self.pix_release(k, deleted);
            }
        }
        self.live_segs -= 1;
        self.dead_segs.push(id);
    }

    /// Live segments within distance 1 of point `p` (a superset).
    fn segs_near_point(&self, p: Point, out: &mut Vec<u32>) {
        out.clear();
        for &t in self.seg_grid.get(cell_of(p, self.sh)) {
            out.push(t);
        }
        out.extend(self.seg_grid.long.iter().copied());
        let segs = &self.segs;
        out.retain(|&t| {
            let s = &segs[t as usize];
            s.live
                && p.x >= s.bb.min.x - 1
                && p.x <= s.bb.max.x + 1
                && p.y >= s.bb.min.y - 1
                && p.y <= s.bb.max.y + 1
        });
    }

    /// Live candidate pixels within distance 1 of segment `s` (a superset).
    fn pix_near_seg(&self, s: u32, out: &mut Vec<u32>) {
        out.clear();
        let sg = &self.segs[s as usize];
        let bb = sg.bb.expand(1);
        let pix = &self.pix;
        self.pix_grid
            .near_seg(sg.a, sg.b, self.sh, &mut |_, items| {
                for &k in items {
                    if bb.contains_point(pix[k as usize].p) {
                        out.push(k);
                    }
                }
            });
    }

    /// Runs the hot-pixel propagation to its fixpoint (least fixpoint above the current
    /// flags). `on_affect(s, state)` is called for every segment that becomes affected.
    fn propagate(
        &mut self,
        pq: &mut Vec<u32>,
        sq: &mut Vec<u32>,
        on_affect: &mut impl FnMut(u32, &mut State),
    ) {
        let mut buf: Vec<u32> = Vec::new();
        loop {
            if let Some(k) = pq.pop() {
                let p = self.pix[k as usize].p;
                self.segs_near_point(p, &mut buf);
                for &t in &buf {
                    let s = &self.segs[t as usize];
                    if !s.affected && rel(s, p) == R::Meets {
                        self.segs[t as usize].affected = true;
                        sq.push(t);
                        on_affect(t, self);
                    }
                }
            } else if let Some(s) = sq.pop() {
                self.pix_near_seg(s, &mut buf);
                let sg = &self.segs[s as usize];
                for &k in &buf {
                    let px = &self.pix[k as usize];
                    if !px.hot && matches!(rel(sg, px.p), R::Own | R::Meets) {
                        self.pix[k as usize].hot = true;
                        pq.push(k);
                    }
                }
            } else {
                break;
            }
        }
    }

    /// The fragment chain of segment `s` (from `a` to `b`), exactly as the batch snap
    /// rounder computes it.
    fn compute_chain(&self, s: u32, sc: &mut Scratch) -> Vec<Point> {
        let sg = &self.segs[s as usize];
        let (a, b) = (sg.a, sg.b);
        self.pix_near_seg(s, &mut sc.ids);
        let mut chain = Vec::new();
        if !sg.affected {
            // Exact position, split at candidate pixels lying on the interior (segment
            // endpoints: crossing pixels on it would have made it affected).
            sc.ord.clear();
            for &k in &sc.ids {
                let p = self.pix[k as usize].p;
                if in_segment_interior(a, b, p) {
                    sc.ord.push((dist2(a, p), p));
                }
            }
            chain.push(a);
            if !sc.ord.is_empty() {
                sc.ord.sort_unstable();
                let mut cur = a;
                for &(_, x) in &sc.ord {
                    if x != cur {
                        chain.push(x);
                        cur = x;
                    }
                }
            }
            chain.push(b);
            return chain;
        }
        // Rerouted through every candidate pixel it meets (all hot), in order.
        let d = sub(b, a);
        sc.ord.clear();
        sc.near.clear();
        for &k in &sc.ids {
            let p = self.pix[k as usize].p;
            match rel(sg, p) {
                R::Meets => sc.ord.push((dot(sub(p, a), d), p)),
                R::Near => sc.near.push(p),
                _ => {}
            }
        }
        sc.ord.sort_unstable();
        sc.ord.dedup();
        let mut poly = Vec::with_capacity(sc.ord.len() + 2);
        poly.push(a);
        poly.extend(sc.ord.iter().map(|x| x.1));
        poly.push(b);
        // Candidate centres lying on a rerouted fragment's interior split it.
        sc.ins.clear();
        for &c in &sc.near {
            for k in 0..poly.len() - 1 {
                if in_segment_interior(poly[k], poly[k + 1], c) {
                    sc.ins.push((k, dist2(poly[k], c), c));
                    break;
                }
            }
        }
        sc.ins.sort_unstable();
        sc.ins.dedup();
        chain.push(a);
        let mut j = 0;
        for k in 0..poly.len() - 1 {
            let mut cur = poly[k];
            while j < sc.ins.len() && sc.ins[j].0 == k {
                if sc.ins[j].2 != cur {
                    chain.push(sc.ins[j].2);
                    cur = sc.ins[j].2;
                }
                j += 1;
            }
            if poly[k + 1] != cur {
                chain.push(poly[k + 1]);
            }
        }
        chain
    }

    // ----- arrangement -------------------------------------------------------------------

    /// Adds (or removes) one fragment `p -> q` of an input segment.
    fn frag(&mut self, p: Point, q: Point, operand: u8, tag: u64, add: bool) {
        let (lo, hi, sign) = if p < q { (p, q, 1i8) } else { (q, p, -1i8) };
        let id = match self.edge_map.get(&(lo, hi)) {
            Some(&id) => id,
            None => {
                if !add {
                    return;
                }
                let rec = Edge {
                    lo,
                    hi,
                    contribs: Vec::new(),
                    delta: [0, 0],
                    tag: 0,
                    below: [0, 0],
                    bdir: 0,
                    ring: NONE,
                    live: true,
                    long: false,
                };
                let id = if let Some(id) = self.edge_free.pop() {
                    self.edges[id as usize] = rec;
                    id
                } else {
                    self.edges.push(rec);
                    self.edges.len() as u32 - 1
                };
                let long = self.edge_grid.add_seg(id, lo, hi, self.sh);
                self.edges[id as usize].long = long;
                self.edge_map.insert((lo, hi), id);
                id
            }
        };
        let e = &mut self.edges[id as usize];
        match e
            .contribs
            .iter()
            .position(|c| c.operand == operand && c.sign == sign && c.tag == tag)
        {
            Some(i) => {
                if add {
                    e.contribs[i].count += 1;
                } else {
                    e.contribs[i].count -= 1;
                    if e.contribs[i].count == 0 {
                        e.contribs.swap_remove(i);
                    }
                }
            }
            None => {
                if !add {
                    return;
                }
                e.contribs.push(Contrib {
                    operand,
                    sign,
                    tag,
                    count: 1,
                });
            }
        }
        // Net winding change per operand, and the tag rule of the batch merge: the first
        // contributor by (operand, sign, tag) whose operand has a non-zero net change.
        let mut delta = [0i32; 2];
        for c in &e.contribs {
            delta[(c.operand & 1) as usize] += c.sign as i32 * c.count as i32;
        }
        e.delta = delta;
        e.tag = e
            .contribs
            .iter()
            .filter(|c| delta[(c.operand & 1) as usize] != 0)
            .map(|c| (c.operand, c.sign, c.tag))
            .min()
            .map_or(0, |x| x.2);
        if e.contribs.is_empty() {
            e.live = false;
            let long = e.long;
            self.edge_map.remove(&(lo, hi));
            self.edge_grid.del_seg(id, lo, hi, self.sh, long);
            self.dead_edges.push(id);
        }
    }

    fn sort_sweep(&self, v: &mut Vec<u32>) {
        let edges = &self.edges;
        let keyed: Vec<(Point, Point, u32)> = v
            .iter()
            .map(|&e| (edges[e as usize].lo, edges[e as usize].hi, e))
            .collect();
        let keyed = crate::par::bucket_sort_by_x(
            keyed,
            |k| k.0.x,
            |a, b| cmp_sweep_edges((a.0, a.1), (b.0, b.1)).then_with(|| a.1.cmp(&b.1)),
        );
        v.clear();
        v.extend(keyed.into_iter().map(|k| k.2));
    }

    /// Sets the winding below edge `e` and derives its boundary direction.
    fn set_winding(&mut self, e: u32, below: [i32; 2]) {
        let d = self.edges[e as usize].delta;
        let above = [below[0] + d[0], below[1] + d[1]];
        let (ib, ia) = (self.inside(below), self.inside(above));
        let ed = &mut self.edges[e as usize];
        ed.below = below;
        ed.bdir = if ib == ia {
            0
        } else if ia {
            1
        } else {
            -1
        };
    }

    /// The edge just below an edge starting at `v` in direction `d` (`d` pointing into the
    /// right half-plane), in the sense of the batch sweep: among the edges starting at `v`
    /// the closest one angularly below `d`, else the topmost edge passing strictly below
    /// `v`. With `bonly`, only boundary edges count. Also returns the rounded-down `y` of
    /// the hit point at `v.x` (`v.y` for an edge at `v`, a large negative value for none).
    fn edge_below(&self, v: Point, d: Point, bonly: bool) -> (u32, i64) {
        let edges = &self.edges;
        let ok = |e: &Edge| e.active() && (!bonly || e.bdir != 0);
        let sh = self.sh;
        let dir = |id: u32| {
            let e = &edges[id as usize];
            sub(e.hi, e.lo)
        };
        let mut best1 = NONE;
        let cell_v = self.edge_grid.get(cell_of(v, sh));
        for &id in cell_v.iter().chain(self.edge_grid.long.iter()) {
            let e = &edges[id as usize];
            if e.lo == v
                && ok(e)
                && cmp_dir_halfplane(dir(id), d) == Ordering::Less
                && (best1 == NONE || cmp_dir_halfplane(dir(best1), dir(id)) == Ordering::Less)
            {
                best1 = id;
            }
        }
        if best1 != NONE {
            return (best1, v.y);
        }
        // Topmost edge strictly below v among those spanning it.
        let better = |id: u32, best: u32| -> u32 {
            let e = &edges[id as usize];
            if !(ok(e) && e.lo < v && v < e.hi && orient(e.lo, e.hi, v) > 0) {
                return best;
            }
            if best == NONE {
                return id;
            }
            let b = &edges[best as usize];
            let (n1, d1) = y_at(e.lo, e.hi, v.x);
            let (n2, d2) = y_at(b.lo, b.hi, v.x);
            // Equal heights: both start at that point; the steeper one is above.
            let c = (n1 * d2)
                .cmp(&(n2 * d1))
                .then_with(|| cmp_dir_halfplane(sub(e.hi, e.lo), sub(b.hi, b.lo)));
            if c == Ordering::Greater { id } else { best }
        };
        let floor_y = |id: u32| -> i64 {
            let e = &edges[id as usize];
            let (n, d) = y_at(e.lo, e.hi, v.x);
            floor_div(n, d) as i64
        };
        let mut best = NONE;
        for &id in &self.edge_grid.long {
            best = better(id, best);
        }
        let cx = v.x >> sh;
        let mut row = v.y >> sh;
        let min_row = self.edge_grid.min_row;
        let ncells = self.edge_grid.occupied as i64;
        while row >= min_row {
            if row - min_row > ncells + 64 {
                // Long empty column: scan the occupied cells of the column instead.
                self.edge_grid.for_each_cell(&mut |c, items| {
                    if c.0 == cx && c.1 <= row {
                        for &id in items {
                            best = better(id, best);
                        }
                    }
                });
                break;
            }
            for &id in self.edge_grid.get((cx, row)) {
                best = better(id, best);
            }
            if best != NONE && (floor_y(best) >> sh) >= row {
                break;
            }
            row -= 1;
        }
        if best == NONE {
            (NONE, -(1i64 << 42))
        } else {
            (best, floor_y(best))
        }
    }

    // ----- rings -------------------------------------------------------------------------

    /// Links boundary edges `bnd` into rings (new ring records, ids appended to `out`).
    /// Fails when the edges are not balanced at every vertex.
    fn link(&mut self, bnd: &[u32], out: &mut Vec<u32>) -> core::result::Result<(), &'static str> {
        let dir: Vec<DirEdge> = bnd
            .iter()
            .map(|&e| {
                let ed = &self.edges[e as usize];
                let (from, to) = ed.ends();
                DirEdge {
                    from,
                    to,
                    tag: ed.tag,
                    below: 0,
                    above: 0,
                }
            })
            .collect();
        // Balance check: in-degree equals out-degree everywhere.
        let mut bal: FxMap<Point, i32> = FxMap::default();
        for d in &dir {
            *bal.entry(d.from).or_insert(0) += 1;
            *bal.entry(d.to).or_insert(0) -= 1;
        }
        if bal.values().any(|&x| x != 0) {
            return Err("unbalanced boundary");
        }
        drop(bal);
        let first_new = out.len();
        let (walks, _) = link_rings(&dir);
        let mut used = 0usize;
        for w in walks {
            used += w.len();
            if w.is_empty() {
                continue;
            }
            let pts: Vec<Point> = w.iter().map(|&k| dir[k as usize].from).collect();
            let tags: Vec<u64> = w.iter().map(|&k| dir[k as usize].tag).collect();
            let area = ring_area2(&pts);
            let real = pts.len() >= 3 && area != 0;
            // Query edge: at the smallest vertex, the lower of the ring's two edges there.
            let m = w.len();
            let (i, _) = pts
                .iter()
                .enumerate()
                .min_by_key(|(_, p)| **p)
                .unwrap_or((0, &Point::default()));
            let v = pts[i];
            let out_e = &dir[w[i] as usize];
            let in_e = &dir[w[(i + m - 1) % m] as usize];
            let d_out = sub(out_e.to, v);
            let d_in = sub(in_e.from, v);
            let qd = if cmp_dir_halfplane(d_out, d_in) == Ordering::Less {
                d_out
            } else {
                d_in
            };
            let mut raw = RawRing { pts, tags };
            if real {
                let mut pinch: Vec<Point> = raw
                    .pts
                    .iter()
                    .copied()
                    .filter(|p| self.outdeg.get(p).is_some_and(|&n| n > 1))
                    .collect();
                pinch.sort_unstable();
                remove_collinear(&mut raw, &pinch);
                rotate_to_min(&mut raw);
            }
            let edges: Vec<u32> = w.iter().map(|&k| bnd[k as usize]).collect();
            let rec = RingRec {
                edges,
                raw,
                is_hole: area < 0,
                real,
                live: true,
                wit: Wit {
                    edge: NONE,
                    v,
                    ylo: v.y,
                },
                query_dir: qd,
                wit_reg: false,
                wit_long: false,
                gnext: NONE,
                stamp: self.next_stamp,
            };
            self.next_stamp += 1;
            let id = if let Some(id) = self.ring_free.pop() {
                self.rings[id as usize] = rec;
                id
            } else {
                self.rings.push(rec);
                self.rings.len() as u32 - 1
            };
            for &e in &self.rings[id as usize].edges {
                self.edges[e as usize].ring = id;
            }
            out.push(id);
        }
        if used != bnd.len() {
            return Err("open walk");
        }
        self.group_rings(&dir, bnd, &out[first_new..]);
        Ok(())
    }

    /// Groups the rings just linked from `dir` (ids `rings`) by the closed walk they were
    /// split from: rings are in the same walk exactly when the linker's successor (first
    /// out-edge clockwise from the reversed incoming edge) leads from one to the other,
    /// which can only happen at vertices with several out-edges.
    fn group_rings(&mut self, dir: &[DirEdge], bnd: &[u32], rings: &[u32]) {
        for &r in rings {
            self.rings[r as usize].gnext = r;
        }
        // Out-edges of vertices with several out-edges, counter-clockwise by direction.
        let d = |k: u32| sub(dir[k as usize].to, dir[k as usize].from);
        let mut cnt: FxMap<Point, u32> = FxMap::default();
        for e in dir {
            if self.outdeg.get(&e.from).is_some_and(|&n| n > 1) {
                *cnt.entry(e.from).or_insert(0) += 1;
            }
        }
        cnt.retain(|_, n| *n > 1);
        if cnt.is_empty() {
            return;
        }
        let mut perm: Vec<u32> = (0..dir.len() as u32)
            .filter(|&k| cnt.contains_key(&dir[k as usize].from))
            .collect();
        perm.sort_unstable_by(|&a, &b| {
            dir[a as usize]
                .from
                .cmp(&dir[b as usize].from)
                .then_with(|| cmp_angle(d(a), d(b)))
        });
        let n = perm.len();
        let mut runs: FxMap<Point, (usize, usize)> = FxMap::default();
        let mut i = 0;
        while i < n {
            let v = dir[perm[i] as usize].from;
            let mut j = i + 1;
            while j < n && dir[perm[j] as usize].from == v {
                j += 1;
            }
            runs.insert(v, (i, j));
            i = j;
        }
        let ring_of = |k: u32| self.edges[bnd[k as usize] as usize].ring;
        let mut uf: FxMap<u32, u32> = FxMap::default();
        fn find(uf: &mut FxMap<u32, u32>, mut x: u32) -> u32 {
            while let Some(&p) = uf.get(&x) {
                if p == x {
                    break;
                }
                let gp = uf.get(&p).copied().unwrap_or(p);
                uf.insert(x, gp);
                x = gp;
            }
            x
        }
        let mut unions: Vec<(u32, u32)> = Vec::new();
        for (k, e) in dir.iter().enumerate() {
            let Some(&(lo, hi)) = runs.get(&e.to) else {
                continue;
            };
            let r = sub(e.from, e.to);
            let at = lo + perm[lo..hi].partition_point(|&k| cmp_angle(d(k), r) == Ordering::Less);
            let next = perm[if at == lo { hi - 1 } else { at - 1 }];
            let (a, b) = (ring_of(k as u32), ring_of(next));
            if a != b && a != NONE && b != NONE {
                unions.push((a, b));
            }
        }
        if unions.is_empty() {
            return;
        }
        for &(a, b) in &unions {
            uf.entry(a).or_insert(a);
            uf.entry(b).or_insert(b);
            let (ra, rb) = (find(&mut uf, a), find(&mut uf, b));
            if ra != rb {
                let (lo, hi) = if ra < rb { (ra, rb) } else { (rb, ra) };
                uf.insert(hi, lo);
            }
        }
        // Circular lists per group.
        let mut members: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
        let keys: Vec<u32> = uf.keys().copied().collect();
        for x in keys {
            let r = find(&mut uf, x);
            members.entry(r).or_default().push(x);
        }
        for (_, m) in members {
            for (i, &r) in m.iter().enumerate() {
                self.rings[r as usize].gnext = m[(i + 1) % m.len()];
            }
        }
    }

    /// Finds the boundary edge below ring `r`'s query edge.
    fn shoot_ring(&mut self, r: u32) {
        self.wit_unregister(r);
        let rr = &self.rings[r as usize];
        let (v, d) = (rr.wit.v, rr.query_dir);
        let (e, ylo) = self.edge_below(v, d, true);
        let long = self.wit_grid.add_seg(r, Point::new(v.x, ylo), v, self.sh);
        let rr = &mut self.rings[r as usize];
        rr.wit = Wit { edge: e, v, ylo };
        rr.wit_reg = true;
        rr.wit_long = long;
        if e != NONE {
            self.wit_of.entry(e).or_default().push(r);
        }
    }

    /// Removes ring `r`'s witness from the indexes.
    fn wit_unregister(&mut self, r: u32) {
        let rr = &mut self.rings[r as usize];
        if !rr.wit_reg {
            return;
        }
        rr.wit_reg = false;
        let (w, long) = (rr.wit, rr.wit_long);
        self.wit_grid
            .del_seg(r, Point::new(w.v.x, w.ylo), w.v, self.sh, long);
        if w.edge != NONE
            && let Some(v) = self.wit_of.get_mut(&w.edge)
        {
            if let Some(i) = v.iter().position(|&x| x == r) {
                v.swap_remove(i);
            }
            if v.is_empty() {
                self.wit_of.remove(&w.edge);
            }
        }
    }

    /// The parent of ring `r` from its witness, given the parents of the rings below.
    /// `Err` when the witness leads back to `r` itself.
    fn parent_from_witness(&self, r: u32) -> core::result::Result<Option<u32>, &'static str> {
        let w = self.rings[r as usize].wit.edge;
        if w == NONE {
            return Ok(None);
        }
        let eb = &self.edges[w as usize];
        let s = eb.ring;
        if s == NONE || !self.rings[s as usize].live || !self.rings[s as usize].real {
            return Ok(None);
        }
        if s == r {
            return Err("nesting cycle");
        }
        let s_hole = self.rings[s as usize].is_hole;
        let need_parent_of_s = if eb.bdir > 0 { s_hole } else { !s_hole };
        Ok(if need_parent_of_s {
            self.parent.get(s as usize).copied().flatten()
        } else {
            Some(s)
        })
    }

    /// Recomputes the parents of the rings in `dirty` and of every ring whose parent
    /// depends on a changed one. `dirty` must contain the new rings and every ring whose
    /// witness edge belongs to a new ring. Rings are processed bottom-up (by their query vertex and
    /// direction), so each ring's dependencies are final when it is evaluated.
    fn update_parents(&mut self, dirty: &[u32]) -> core::result::Result<(), &'static str> {
        #[derive(PartialEq, Eq)]
        struct Key(Point, Point, u32);
        impl Ord for Key {
            fn cmp(&self, o: &Self) -> Ordering {
                self.0
                    .cmp(&o.0)
                    .then_with(|| cmp_dir_halfplane(self.1, o.1))
                    .then_with(|| self.2.cmp(&o.2))
            }
        }
        impl PartialOrd for Key {
            fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
                Some(self.cmp(o))
            }
        }
        if self.parent.len() < self.rings.len() {
            self.parent.resize(self.rings.len(), None);
        }
        let ep = self.next_epoch();
        let mut queue: std::collections::BTreeSet<Key> = std::collections::BTreeSet::new();
        let key = |st: &State, r: u32| {
            Key(
                st.rings[r as usize].wit.v,
                st.rings[r as usize].query_dir,
                r,
            )
        };
        for &r in dirty {
            let rr = &self.rings[r as usize];
            if rr.live && rr.real && self.ring_mark[r as usize] != ep {
                self.ring_mark[r as usize] = ep;
                queue.insert(key(self, r));
            }
        }
        let mut last: Option<Key> = None;
        while let Some(k) = queue.pop_first() {
            if last.as_ref().is_some_and(|l| *l >= k) {
                return Err("nesting order");
            }
            let r = k.2;
            self.ring_mark[r as usize] = 0;
            let p = self.parent_from_witness(r)?;
            let changed = self.parent[r as usize] != p;
            self.parent[r as usize] = p;
            last = Some(k);
            // Dependents of a new ring (whose stale parent slot means nothing) are in
            // `dirty` already; others only need a look when this parent changed.
            if !changed {
                continue;
            }
            // Rings whose parent is inherited from r's.
            let is_hole = self.rings[r as usize].is_hole;
            for i in 0..self.rings[r as usize].edges.len() {
                let e = self.rings[r as usize].edges[i];
                let inherit = if self.edges[e as usize].bdir > 0 {
                    is_hole
                } else {
                    !is_hole
                };
                if !inherit {
                    continue;
                }
                let Some(deps) = self.wit_of.get(&e) else {
                    continue;
                };
                for q in deps.clone() {
                    let qr = &self.rings[q as usize];
                    if qr.live && qr.real && self.ring_mark[q as usize] != ep {
                        self.ring_mark[q as usize] = ep;
                        queue.insert(key(self, q));
                    }
                }
            }
        }
        Ok(())
    }

    /// Parents of all rings, from their witnesses (the batch nesting rule).
    fn compute_parents(&mut self) -> core::result::Result<(), &'static str> {
        let n = self.rings.len();
        self.parent.clear();
        self.parent.resize(n, None);
        let mut state = vec![0u8; n];
        let mut stack: Vec<u32> = Vec::new();
        for r0 in 0..n {
            if !self.rings[r0].live || !self.rings[r0].real || state[r0] != 0 {
                continue;
            }
            state[r0] = 1;
            stack.push(r0 as u32);
            while let Some(&r) = stack.last() {
                let rr = &self.rings[r as usize];
                let w = rr.wit.edge;
                let res: Option<Option<u32>> = if w == NONE {
                    Some(None)
                } else {
                    let eb = &self.edges[w as usize];
                    let s = eb.ring;
                    if s == NONE || !self.rings[s as usize].live || !self.rings[s as usize].real {
                        Some(None)
                    } else {
                        let s_hole = self.rings[s as usize].is_hole;
                        let need_parent_of_s = if eb.bdir > 0 { s_hole } else { !s_hole };
                        if !need_parent_of_s {
                            Some(Some(s))
                        } else if state[s as usize] == 2 {
                            Some(self.parent[s as usize])
                        } else if state[s as usize] == 1 {
                            return Err("nesting cycle");
                        } else {
                            state[s as usize] = 1;
                            stack.push(s);
                            None
                        }
                    }
                };
                if let Some(p) = res {
                    self.parent[r as usize] = p;
                    state[r as usize] = 2;
                    stack.pop();
                }
            }
        }
        Ok(())
    }

    /// The canonical tree of the current rings.
    /// Every outer ring with its holes (sorted), polygons sorted by outer ring: the
    /// canonical flattening of the tree, as [`PolyTree::to_polygon_set`] does it.
    fn flatten(&self) -> Vec<(&RawRing, Vec<&RawRing>)> {
        let n = self.rings.len();
        let ok = |i: usize| self.rings[i].live && self.rings[i].real;
        let mut holes: Vec<Vec<u32>> = vec![Vec::new(); n];
        let mut outers: Vec<u32> = Vec::new();
        for i in 0..n {
            if !ok(i) {
                continue;
            }
            if !self.rings[i].is_hole {
                outers.push(i as u32);
            } else if let Some(p) = self.parent.get(i).copied().flatten()
                && ok(p as usize)
            {
                holes[p as usize].push(i as u32);
            }
        }
        let pts = |i: &u32| &self.rings[*i as usize].raw.pts;
        outers.sort_by(|a, b| pts(a).cmp(pts(b)));
        outers
            .iter()
            .map(|&o| {
                let mut h = core::mem::take(&mut holes[o as usize]);
                h.sort_by(|a, b| pts(a).cmp(pts(b)));
                (
                    &self.rings[o as usize].raw,
                    h.iter().map(|&k| &self.rings[k as usize].raw).collect(),
                )
            })
            .collect()
    }

    /// The canonical tree of the live rings (as `canonical_tree` builds it), with the ring
    /// slot and stamp of every node. Rings found in `old` (same slot and stamp) are moved
    /// from there rather than copied, and its order presorts the children lists.
    fn tree(&self, old: Option<(PolyTree, Vec<(u32, u64)>)>) -> (PolyTree, Vec<(u32, u64)>) {
        let (mut old_nodes, old_meta) = match old {
            Some((t, meta)) if meta.len() == t.nodes.len() => (t.nodes, meta),
            _ => (Vec::new(), Vec::new()),
        };
        let slots = if old_meta.is_empty() {
            0
        } else {
            self.rings.len()
        };
        let mut node_of_slot = vec![NONE; slots];
        for (n, &(slot, _)) in old_meta.iter().enumerate() {
            if let Some(x) = node_of_slot.get_mut(slot as usize) {
                *x = n as u32;
            }
        }
        let mut idx = vec![NONE; self.rings.len()];
        let mut rings: Vec<RawRing> = Vec::new();
        let mut is_hole: Vec<bool> = Vec::new();
        let mut ids: Vec<u32> = Vec::new();
        // Position of every ring in the old tree (`NONE` for new rings).
        let mut hint: Vec<u32> = Vec::new();
        for (i, r) in self.rings.iter().enumerate() {
            if r.live && r.real {
                idx[i] = rings.len() as u32;
                let n = node_of_slot.get(i).copied().unwrap_or(NONE);
                if n != NONE && old_meta[n as usize].1 == r.stamp {
                    let node = &mut old_nodes[n as usize];
                    rings.push(RawRing {
                        pts: core::mem::take(&mut node.ring.0),
                        tags: core::mem::take(&mut node.tags),
                    });
                    hint.push(n);
                } else {
                    rings.push(r.raw.clone());
                    hint.push(NONE);
                }
                is_hole.push(r.is_hole);
                ids.push(i as u32);
            }
        }
        drop(old_nodes);
        let parent: Vec<Option<u32>> = ids
            .iter()
            .map(|&i| {
                self.parent
                    .get(i as usize)
                    .copied()
                    .flatten()
                    .and_then(|p| (idx[p as usize] != NONE).then_some(idx[p as usize]))
            })
            .collect();
        // `canonical_tree`: children (and roots) sorted stably by vertex sequence, that is
        // by (vertex sequence, index); then depth first. Presorting by the old positions
        // leaves the (comparison-heavy) sort little to do.
        let m = rings.len();
        let mut children: Vec<Vec<u32>> = vec![Vec::new(); m];
        let mut roots: Vec<u32> = Vec::new();
        for (r, par) in parent.iter().enumerate() {
            match par {
                Some(p) => children[*p as usize].push(r as u32),
                None => roots.push(r as u32),
            }
        }
        let presort = !old_meta.is_empty();
        let order_list = |l: &mut Vec<u32>| {
            if l.len() < 2 {
                return;
            }
            if presort {
                l.sort_unstable_by_key(|&r| hint[r as usize]);
            }
            l.sort_by(|&a, &b| {
                rings[a as usize]
                    .pts
                    .cmp(&rings[b as usize].pts)
                    .then(a.cmp(&b))
            });
        };
        order_list(&mut roots);
        for c in children.iter_mut() {
            order_list(c);
        }
        let mut new_id = vec![0usize; m];
        let mut order: Vec<u32> = Vec::with_capacity(m);
        let mut stack: Vec<u32> = roots.iter().rev().copied().collect();
        while let Some(r) = stack.pop() {
            new_id[r as usize] = order.len();
            order.push(r);
            stack.extend(children[r as usize].iter().rev());
        }
        let mut taken: Vec<Option<RawRing>> = rings.into_iter().map(Some).collect();
        let nodes = order
            .iter()
            .map(|&r| {
                let rr = taken[r as usize].take().unwrap_or(RawRing {
                    pts: Vec::new(),
                    tags: Vec::new(),
                });
                PolyNode {
                    ring: Ring(rr.pts),
                    tags: rr.tags,
                    is_hole: is_hole[r as usize],
                    parent: parent[r as usize].map(|p| new_id[p as usize]),
                    children: children[r as usize]
                        .iter()
                        .map(|&c| new_id[c as usize])
                        .collect(),
                }
            })
            .collect();
        let meta = order
            .iter()
            .map(|&r| {
                let slot = ids[r as usize];
                (slot, self.rings[slot as usize].stamp)
            })
            .collect();
        (
            PolyTree {
                nodes,
                roots: roots.iter().map(|&r| new_id[r as usize]).collect(),
            },
            meta,
        )
    }

    // ----- incremental update ------------------------------------------------------------

    /// Removes segments `rem`, inserts `add`; `rects` are the bounding boxes of the removed
    /// and inserted rings. Returns the ids of the inserted segments, or `Err` when an
    /// internal consistency check failed (the caller then rebuilds).
    fn apply(
        &mut self,
        rem: &[u32],
        add: &[NewSeg],
        mut d: Vec<Rect>,
    ) -> core::result::Result<Vec<u32>, &'static str> {
        let sh = self.sh;
        let mut buf: Vec<u32> = Vec::new();
        let mut pts: Vec<Point> = Vec::new();

        // 1. Rounding clusters touched by the removed segments (old state): flood the
        //    realized "affected segment meets hot pixel" graph from them.
        let ep = self.next_epoch();
        let mut x_segs: Vec<u32> = Vec::new();
        let mut x_pix: Vec<u32> = Vec::new();
        let mut fq_s: Vec<u32> = Vec::new();
        let mut fq_p: Vec<u32> = Vec::new();
        for &r in rem {
            let s = &self.segs[r as usize];
            if !s.live {
                continue;
            }
            if s.affected && self.seg_mark[r as usize] != ep {
                self.seg_mark[r as usize] = ep;
                fq_s.push(r);
            }
            self.crossings(r, &mut pts);
            pts.push(s.a);
            pts.push(s.b);
            for &p in &pts {
                if let Some(&k) = self.pix_map.get(&p)
                    && self.pix[k as usize].hot
                    && self.pix_mark[k as usize] != ep
                {
                    self.pix_mark[k as usize] = ep;
                    fq_p.push(k);
                }
            }
        }
        loop {
            if let Some(s) = fq_s.pop() {
                x_segs.push(s);
                self.pix_near_seg(s, &mut buf);
                let sg = &self.segs[s as usize];
                for &k in &buf {
                    let px = &self.pix[k as usize];
                    if px.hot
                        && self.pix_mark[k as usize] != ep
                        && matches!(rel(sg, px.p), R::Own | R::Meets)
                    {
                        self.pix_mark[k as usize] = ep;
                        fq_p.push(k);
                    }
                }
            } else if let Some(k) = fq_p.pop() {
                x_pix.push(k);
                let p = self.pix[k as usize].p;
                self.segs_near_point(p, &mut buf);
                for &t in &buf {
                    let s = &self.segs[t as usize];
                    if s.affected
                        && self.seg_mark[t as usize] != ep
                        && matches!(rel(s, p), R::Own | R::Meets)
                    {
                        self.seg_mark[t as usize] = ep;
                        fq_s.push(t);
                    }
                }
            } else {
                break;
            }
        }

        // 2. Remove and insert segments (candidate pixel counts).
        let mut p_ch: Vec<Point> = Vec::new();
        let mut rem_chains: Vec<(u8, u64, Vec<Point>)> = Vec::with_capacity(rem.len());
        for &r in rem {
            if !self.segs[r as usize].live {
                continue;
            }
            self.unregister_seg(r, &mut p_ch);
            let s = &mut self.segs[r as usize];
            rem_chains.push((s.operand, s.tag, core::mem::take(&mut s.chain)));
        }
        let mut new_cross: Vec<u32> = Vec::new();
        let mut new_ids: Vec<u32> = Vec::with_capacity(add.len());
        for s in add {
            let id = self.alloc_seg(s);
            self.register_seg(id, &mut p_ch, &mut new_cross);
            new_ids.push(id);
        }

        // 3. Reset the touched clusters and re-run the propagation.
        let ep2 = self.next_epoch();
        // Segments whose affected flag may change, with their old flag.
        let mut aff_log: Vec<(u32, bool)> = Vec::new();
        for &s in &x_segs {
            if self.segs[s as usize].live {
                self.segs[s as usize].affected = false;
                self.seg_mark[s as usize] = ep2;
                aff_log.push((s, true));
            }
        }
        for &k in &x_pix {
            if self.pix[k as usize].live {
                self.pix[k as usize].hot = false;
            }
        }
        for &s in &new_ids {
            self.seg_mark[s as usize] = ep2;
            aff_log.push((s, false));
        }
        let mut pq: Vec<u32> = Vec::new();
        let mut sq: Vec<u32> = Vec::new();
        let activate = |st: &mut State, k: u32, pq: &mut Vec<u32>| {
            let px = &mut st.pix[k as usize];
            if px.live && !px.hot {
                px.hot = true;
                pq.push(k);
            }
        };
        for &k in x_pix.iter().chain(new_cross.iter()) {
            if self.pix[k as usize].live && self.pix[k as usize].cross > 0 {
                activate(self, k, &mut pq);
            }
        }
        // New candidate pixels met by affected segments.
        for &p in &p_ch {
            let Some(&k) = self.pix_map.get(&p) else {
                continue;
            };
            self.segs_near_point(p, &mut buf);
            for &t in &buf {
                let s = &self.segs[t as usize];
                if s.affected && matches!(rel(s, p), R::Own | R::Meets) {
                    activate(self, k, &mut pq);
                    break;
                }
            }
        }
        // New segments meeting hot pixels.
        for &s in &new_ids {
            self.pix_near_seg(s, &mut buf);
            let sg = &self.segs[s as usize];
            if buf.iter().any(|&k| {
                let px = &self.pix[k as usize];
                px.hot && rel(sg, px.p) == R::Meets
            }) {
                self.segs[s as usize].affected = true;
                sq.push(s);
            }
        }
        let mut log2: Vec<u32> = Vec::new();
        self.propagate(&mut pq, &mut sq, &mut |t, st: &mut State| {
            if st.seg_mark[t as usize] != ep2 {
                st.seg_mark[t as usize] = ep2;
                log2.push(t);
            }
        });
        aff_log.extend(log2.into_iter().map(|t| (t, false)));

        // 4. Segments whose fragments may change.
        let ep3 = self.next_epoch();
        let mut dirty: Vec<u32> = Vec::new();
        for &s in &new_ids {
            self.seg_mark[s as usize] = ep3;
        }
        for &(s, old) in &aff_log {
            if self.segs[s as usize].live
                && self.segs[s as usize].affected != old
                && self.seg_mark[s as usize] != ep3
            {
                self.seg_mark[s as usize] = ep3;
                dirty.push(s);
            }
        }
        for &p in &p_ch {
            self.segs_near_point(p, &mut buf);
            for &t in &buf {
                if self.seg_mark[t as usize] != ep3 {
                    self.seg_mark[t as usize] = ep3;
                    dirty.push(t);
                }
            }
        }
        let mut sc = Scratch::default();
        let mut changed: Vec<(u32, Vec<Point>)> = Vec::new();
        for &s in &dirty {
            let chain = self.compute_chain(s, &mut sc);
            if chain != self.segs[s as usize].chain {
                if let Some(r) = chain_diff_rect(&self.segs[s as usize].chain, &chain) {
                    d.push(r);
                }
                changed.push((s, chain));
            }
        }
        let new_chains: Vec<Vec<Point>> = new_ids
            .iter()
            .map(|&s| self.compute_chain(s, &mut sc))
            .collect();
        if d.is_empty() {
            for (&s, c) in new_ids.iter().zip(new_chains) {
                self.segs[s as usize].chain = c;
            }
            self.finish_update();
            return Ok(new_ids);
        }

        // 5. Old arrangement edges meeting D.
        let ep4 = self.next_epoch();
        let mut old_d: Vec<u32> = Vec::new();
        self.edges_meeting(&d, ep4, false, &mut old_d);
        let old_b: Vec<(u32, Point, Point, u32)> = old_d
            .iter()
            .filter(|&&e| self.edges[e as usize].bdir != 0)
            .map(|&e| {
                let ed = &self.edges[e as usize];
                let (f, t) = ed.ends();
                (e, f, t, ed.ring)
            })
            .collect();

        // 6. Update the fragments.
        for (op, tag, chain) in &rem_chains {
            for w in chain.windows(2) {
                self.frag(w[0], w[1], *op, *tag, false);
            }
        }
        for (s, chain) in changed {
            let sg = &self.segs[s as usize];
            let (op, tag) = (sg.operand, sg.tag);
            let old = core::mem::take(&mut self.segs[s as usize].chain);
            // Only the differing fragments: the others may lie outside D and must keep
            // their records (and windings).
            let (gone, came) = chain_diff(&old, &chain);
            for (p, q) in gone {
                self.frag(p, q, op, tag, false);
            }
            for (p, q) in came {
                self.frag(p, q, op, tag, true);
            }
            self.segs[s as usize].chain = chain;
        }
        for (&s, chain) in new_ids.iter().zip(new_chains) {
            let sg = &self.segs[s as usize];
            let (op, tag) = (sg.operand, sg.tag);
            for w in chain.windows(2) {
                self.frag(w[0], w[1], op, tag, true);
            }
            self.segs[s as usize].chain = chain;
        }
        self.edge_mark.resize(self.edges.len(), 0);
        self.edge_mark2.resize(self.edges.len(), 0);

        // 7. New edges meeting D: windings in sweep order.
        let mut new_d: Vec<u32> = Vec::new();
        self.edges_meeting(&d, ep4, true, &mut new_d);
        self.sort_sweep(&mut new_d);
        for &e in &old_d {
            let ed = &mut self.edges[e as usize];
            ed.bdir = 0;
            ed.ring = NONE;
        }
        for &e in &new_d {
            let ed = &self.edges[e as usize];
            let (lo, dd) = (ed.lo, sub(ed.hi, ed.lo));
            let (f, _) = self.edge_below(lo, dd, false);
            let below = if f == NONE {
                [0, 0]
            } else {
                let fe = &self.edges[f as usize];
                [fe.below[0] + fe.delta[0], fe.below[1] + fe.delta[1]]
            };
            self.set_winding(e, below);
            self.edges[e as usize].ring = NONE;
        }
        let new_b: Vec<u32> = new_d
            .iter()
            .copied()
            .filter(|&e| self.edges[e as usize].bdir != 0)
            .collect();

        // 8. Pinch degrees; dirty vertices; touched rings.
        for &(_, f, _, _) in &old_b {
            if let Some(n) = self.outdeg.get_mut(&f) {
                *n -= 1;
                if *n == 0 {
                    self.outdeg.remove(&f);
                }
            }
        }
        for &e in &new_b {
            let (f, _) = self.edges[e as usize].ends();
            *self.outdeg.entry(f).or_insert(0) += 1;
        }
        let ep5 = self.next_epoch();
        let mut touched: Vec<u32> = Vec::new();
        let touch = |st: &mut State, r: u32, touched: &mut Vec<u32>| {
            let mut x = r;
            while x != NONE && st.ring_mark[x as usize] != ep5 {
                st.ring_mark[x as usize] = ep5;
                touched.push(x);
                x = st.rings[x as usize].gnext;
            }
        };
        for &(_, _, _, r) in &old_b {
            touch(self, r, &mut touched);
        }
        let mut dv: Vec<Point> = Vec::with_capacity(2 * (old_b.len() + new_b.len()));
        for &(_, f, t, _) in &old_b {
            dv.push(f);
            dv.push(t);
        }
        for &e in &new_b {
            let ed = &self.edges[e as usize];
            dv.push(ed.lo);
            dv.push(ed.hi);
        }
        dv.sort_unstable();
        dv.dedup();
        for &v in &dv {
            let c = cell_of(v, sh);
            let n_cell = self.edge_grid.get(c).len();
            for i in 0..n_cell + self.edge_grid.long.len() {
                let e = if i < n_cell {
                    self.edge_grid.get(c)[i]
                } else {
                    self.edge_grid.long[i - n_cell]
                };
                let ed = &self.edges[e as usize];
                if ed.active() && ed.bdir != 0 && (ed.lo == v || ed.hi == v) {
                    let r = ed.ring;
                    touch(self, r, &mut touched);
                }
            }
        }

        // 9. Relink: unchanged edges of touched rings plus new boundary edges.
        let mut free: Vec<u32> = Vec::new();
        for &r in &touched {
            let rr = &self.rings[r as usize];
            for &e in &rr.edges {
                // Edges meeting D are replaced by the new boundary edges below.
                if self.edge_mark[e as usize] != ep4 && self.edge_mark2[e as usize] != ep4 {
                    free.push(e);
                }
            }
        }
        free.extend_from_slice(&new_b);
        // The batch assembler links edges in sweep order, so every closed walk starts at its
        // first edge in that order; splitting a walk at repeated (pinch) vertices depends on
        // where it starts, so keep the same order.
        self.sort_sweep(&mut free);
        for &r in &touched {
            self.wit_unregister(r);
            let rr = &mut self.rings[r as usize];
            rr.live = false;
            rr.edges.clear();
            self.ring_free.push(r);
        }
        let mut new_rings: Vec<u32> = Vec::new();
        self.link(&free, &mut new_rings)?;

        // 10. Nesting: re-shoot rings whose ray meets a changed boundary edge.
        let mut ch: Vec<(Point, Point)> = old_b.iter().map(|&(_, f, t, _)| (f, t)).collect();
        for &e in &new_b {
            let ed = &self.edges[e as usize];
            ch.push((ed.lo, ed.hi));
        }
        let ep6 = self.next_epoch();
        for &r in &new_rings {
            self.ring_mark[r as usize] = ep6;
        }
        let mut dirty_rings: Vec<u32> = new_rings.clone();
        let mut cand: Vec<u32> = Vec::new();
        for &(a, b) in &ch {
            cand.clear();
            self.wit_grid
                .near_seg(a, b, sh, &mut |_, items| cand.extend_from_slice(items));
            cand.extend_from_slice(&self.wit_grid.long);
            for &r in &cand {
                let rr = &self.rings[r as usize];
                if !rr.live || self.ring_mark[r as usize] == ep6 {
                    continue;
                }
                let w = rr.wit;
                if segments_intersect(Point::new(w.v.x, w.ylo), w.v, a, b) {
                    self.ring_mark[r as usize] = ep6;
                    dirty_rings.push(r);
                }
            }
        }
        for &r in &dirty_rings {
            self.shoot_ring(r);
        }
        // Rings whose witness edge now belongs to a new ring.
        for &r in &new_rings {
            for i in 0..self.rings[r as usize].edges.len() {
                let e = self.rings[r as usize].edges[i];
                if let Some(deps) = self.wit_of.get(&e) {
                    dirty_rings.extend_from_slice(deps);
                }
            }
        }
        self.update_parents(&dirty_rings)?;
        self.finish_update();
        Ok(new_ids)
    }

    /// Active edges meeting any rectangle of `d`, each once (marked with `ep` in
    /// `edge_mark2` when `second`, else `edge_mark`).
    fn edges_meeting(&mut self, d: &[Rect], ep: u32, second: bool, out: &mut Vec<u32>) {
        let sh = self.sh;
        let edges = &self.edges;
        let mark = if second {
            &mut self.edge_mark2
        } else {
            &mut self.edge_mark
        };
        mark.resize(edges.len(), 0);
        let mut visit = |id: u32, r: &Rect, out: &mut Vec<u32>| {
            let e = &edges[id as usize];
            if mark[id as usize] != ep && e.active() && segment_meets_rect(e.lo, e.hi, r) {
                mark[id as usize] = ep;
                out.push(id);
            }
        };
        for r in d {
            self.edge_grid.in_rect(r, sh, &mut |items| {
                for &id in items {
                    visit(id, r, out);
                }
            });
            for &id in &self.edge_grid.long {
                visit(id, r, out);
            }
        }
    }

    /// Compares rings and nesting with a fresh build of the same segments.
    fn check_rings(&self) -> Option<String> {
        let input: Vec<NewSeg> = self
            .segs
            .iter()
            .filter(|s| s.live)
            .map(|s| NewSeg {
                a: s.a,
                b: s.b,
                tag: s.tag,
                operand: s.operand,
            })
            .collect();
        let mut fresh = State::empty(self.rule, 0, None);
        State::build(self.rule, &input, &mut fresh);
        let key = |st: &State| {
            let mut v: Vec<(Vec<Point>, bool, Option<Vec<Point>>)> = st
                .rings
                .iter()
                .enumerate()
                .filter(|(_, r)| r.live && r.real)
                .map(|(i, r)| {
                    (
                        r.raw.pts.clone(),
                        r.is_hole,
                        st.parent
                            .get(i)
                            .copied()
                            .flatten()
                            .map(|p| st.rings[p as usize].raw.pts.clone()),
                    )
                })
                .collect();
            v.sort();
            v
        };
        let (a, b) = (key(self), key(&fresh));
        if a != b {
            let only_a: Vec<_> = a.iter().filter(|x| !b.contains(x)).collect();
            let only_b: Vec<_> = b.iter().filter(|x| !a.contains(x)).collect();
            return Some(format!(
                "rings differ:\n  incremental only: {only_a:?}\n  fresh only: {only_b:?}"
            ));
        }
        None
    }

    /// Recomputes the hot-pixel fixpoint and all chains from scratch and compares.
    fn check_snap(&self) -> Option<String> {
        let mut st = self.clone();
        for p in st.pix.iter_mut() {
            p.hot = false;
        }
        for s in st.segs.iter_mut() {
            s.affected = false;
        }
        let mut pq = Vec::new();
        let mut sq = Vec::new();
        for k in 0..st.pix.len() {
            if st.pix[k].live && st.pix[k].cross > 0 {
                st.pix[k].hot = true;
                pq.push(k as u32);
            }
        }
        st.propagate(&mut pq, &mut sq, &mut |_, _| {});
        for k in 0..st.pix.len() {
            if st.pix[k].live && st.pix[k].hot != self.pix[k].hot {
                return Some(format!(
                    "pixel {:?} hot {} want {}",
                    st.pix[k].p, self.pix[k].hot, st.pix[k].hot
                ));
            }
        }
        let mut sc = Scratch::default();
        for i in 0..st.segs.len() {
            if !st.segs[i].live {
                continue;
            }
            if st.segs[i].affected != self.segs[i].affected {
                return Some(format!(
                    "seg {i} {:?}-{:?} affected {} want {}",
                    st.segs[i].a, st.segs[i].b, self.segs[i].affected, st.segs[i].affected
                ));
            }
            let c = st.compute_chain(i as u32, &mut sc);
            if c != self.segs[i].chain {
                return Some(format!(
                    "seg {i} chain {:?} want {:?}",
                    self.segs[i].chain, c
                ));
            }
        }
        None
    }

    /// Recomputes all windings with one sweep and reports the first mismatch.
    fn check_windings(&self) -> Option<String> {
        let mut order: Vec<u32> = (0..self.edges.len() as u32)
            .filter(|&e| self.edges[e as usize].active())
            .collect();
        self.sort_sweep(&mut order);
        let segs: Vec<(Point, Point)> = order
            .iter()
            .map(|&e| (self.edges[e as usize].lo, self.edges[e as usize].hi))
            .collect();
        let mut below = vec![[0i32; 2]; order.len()];
        let mut bel_id = vec![NONE; order.len()];
        sweep(&segs, |k, b| {
            if let Some(b) = b {
                let eb = &self.edges[order[b as usize] as usize];
                below[k as usize] = [
                    below[b as usize][0] + eb.delta[0],
                    below[b as usize][1] + eb.delta[1],
                ];
                bel_id[k as usize] = order[b as usize];
            }
        });
        for (k, &e) in order.iter().enumerate() {
            let ed = &self.edges[e as usize];
            if ed.below != below[k] {
                let (f, _) = self.edge_below(ed.lo, sub(ed.hi, ed.lo), false);
                return Some(format!(
                    "edge {e} {:?}-{:?} below {:?} want {:?}; sweep below {:?}, shoot {:?}",
                    ed.lo,
                    ed.hi,
                    ed.below,
                    below[k],
                    (bel_id[k] != NONE).then(|| {
                        let b = &self.edges[bel_id[k] as usize];
                        (bel_id[k], b.lo, b.hi)
                    }),
                    (f != NONE).then(|| {
                        let b = &self.edges[f as usize];
                        (f, b.lo, b.hi)
                    }),
                ));
            }
        }
        None
    }

    fn finish_update(&mut self) {
        for k in core::mem::take(&mut self.dead_pix) {
            self.pix_free.push(k);
        }
        for e in core::mem::take(&mut self.dead_edges) {
            if !self.edges[e as usize].live {
                self.edge_free.push(e);
            }
        }
        for s in core::mem::take(&mut self.dead_segs) {
            self.seg_free.push(s);
        }
    }
}

/// Scratch buffers for fragment computation.
#[derive(Default)]
struct Scratch {
    ids: Vec<u32>,
    ord: Vec<(i128, Point)>,
    near: Vec<Point>,
    ins: Vec<(usize, i128, Point)>,
}

/// Fragments present only in `old` and only in `new` (as multisets).
#[allow(clippy::type_complexity)]
fn chain_diff(old: &[Point], new: &[Point]) -> (Vec<(Point, Point)>, Vec<(Point, Point)>) {
    let mut a: Vec<(Point, Point)> = old.windows(2).map(|w| (w[0], w[1])).collect();
    let mut b: Vec<(Point, Point)> = new.windows(2).map(|w| (w[0], w[1])).collect();
    a.sort_unstable();
    b.sort_unstable();
    let (mut ga, mut gb) = (Vec::new(), Vec::new());
    let (mut i, mut j) = (0, 0);
    while i < a.len() || j < b.len() {
        match (a.get(i), b.get(j)) {
            (Some(x), Some(y)) if x == y => {
                i += 1;
                j += 1;
            }
            (Some(x), Some(y)) if x < y => {
                ga.push(*x);
                i += 1;
            }
            (Some(_), Some(y)) => {
                gb.push(*y);
                j += 1;
            }
            (Some(x), None) => {
                ga.push(*x);
                i += 1;
            }
            (None, Some(y)) => {
                gb.push(*y);
                j += 1;
            }
            (None, None) => break,
        }
    }
    (ga, gb)
}

/// Bounding box of the fragments present in exactly one of two chains (as multisets).
fn chain_diff_rect(old: &[Point], new: &[Point]) -> Option<Rect> {
    let (ga, gb) = chain_diff(old, new);
    ga.iter()
        .chain(gb.iter())
        .map(|&(p, q)| Rect::new(p, q))
        .reduce(|a, b| a.union(&b))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boolean::{Boolean, Op, boolean};
    use crate::geom::MAX_COORD;

    fn sq(x: i64, y: i64, s: i64) -> Ring {
        Ring::from([(x, y), (x + s, y), (x + s, y + s), (x, y + s)])
    }

    fn reference(zone: &Ring, obst: &[Ring]) -> PolygonSet {
        boolean(Op::Difference, zone, obst, FillRule::NonZero).unwrap()
    }

    fn reference_tree(zone: &Ring, obst: &[Ring]) -> PolyTree {
        Boolean::new()
            .subject(zone, FillRule::NonZero)
            .clip(obst, FillRule::NonZero)
            .op(Op::Difference)
            .execute_tree()
            .unwrap()
    }

    #[test]
    fn empty_engine_is_the_zone() {
        let zone = sq(0, 0, 100);
        let mut z = ZoneFill::new(&zone, FillRule::NonZero).unwrap();
        assert!(z.is_empty());
        assert_eq!(z.fill(), reference(&zone, &[]));
        assert_eq!(z.fill_rule(), FillRule::NonZero);
        z.verify().unwrap();
        let mut e = ZoneFill::new(&Vec::<Ring>::new(), FillRule::NonZero).unwrap();
        e.insert(1, &sq(0, 0, 5)).unwrap();
        assert!(e.fill().is_empty());
    }

    #[test]
    fn bookkeeping() {
        let zone = sq(0, 0, 100);
        let mut z = ZoneFill::new(&zone, FillRule::NonZero).unwrap();
        assert!(!z.insert(3, &sq(10, 10, 5)).unwrap());
        assert!(z.insert(3, &sq(20, 10, 5)).unwrap());
        assert!(!z.update(4, &sq(30, 10, 5)).unwrap());
        assert!(!z.contains(4));
        assert!(z.insert(1, &sq(40, 10, 5)).is_ok());
        assert_eq!(z.ids(), vec![1, 3]);
        assert_eq!(z.len(), 2);
        assert_eq!(z.fill(), reference(&zone, &[sq(40, 10, 5), sq(20, 10, 5)]));
        assert!(z.update(1, &sq(40, 40, 5)).unwrap());
        assert!(z.remove(3));
        assert!(!z.remove(3));
        assert_eq!(z.ids(), vec![1]);
        assert_eq!(z.fill(), reference(&zone, &[sq(40, 40, 5)]));
        // Insert then remove before a commit: nothing happens.
        z.insert(9, &sq(1, 1, 3)).unwrap();
        assert!(z.remove(9));
        assert_eq!(z.len(), 1);
        z.clear();
        assert!(z.is_empty());
        assert_eq!(z.fill(), reference(&zone, &[]));
        z.verify().unwrap();
    }

    #[test]
    fn errors_leave_the_engine_unchanged() {
        let zone = sq(0, 0, 100);
        let bad = Ring::from([(0, 0), (MAX_COORD + 1, 0), (0, 5)]);
        assert!(matches!(
            ZoneFill::new(&bad, FillRule::NonZero),
            Err(Error::CoordinateOutOfRange(_))
        ));
        let mut z = ZoneFill::new(&zone, FillRule::NonZero).unwrap();
        z.insert(1, &sq(10, 10, 10)).unwrap();
        assert!(z.insert(2, &bad).is_err());
        assert!(z.update(1, &bad).is_err());
        assert!(z.set_zone(&bad).is_err());
        assert_eq!(z.ids(), vec![1]);
        assert_eq!(z.fill(), reference(&zone, &[sq(10, 10, 10)]));
    }

    #[test]
    fn identical_reinsert_is_free() {
        let zone = sq(0, 0, 100);
        let mut z = ZoneFill::new(&zone, FillRule::NonZero).unwrap();
        z.set_auto_rebuild(false);
        z.insert(1, &sq(10, 10, 10)).unwrap();
        let t1 = z.fill_tree();
        z.update(1, &sq(10, 10, 10)).unwrap();
        assert_eq!(z.fill_tree(), t1);
        assert_eq!(z.rebuild_count(), 1);
    }

    #[test]
    fn set_zone_keeps_obstacles() {
        let mut z = ZoneFill::new(&sq(0, 0, 100), FillRule::NonZero).unwrap();
        z.insert(1, &sq(10, 10, 10)).unwrap();
        z.insert(2, &sq(150, 10, 10)).unwrap();
        z.set_zone(&sq(0, 0, 200)).unwrap();
        assert_eq!(
            z.fill(),
            reference(&sq(0, 0, 200), &[sq(10, 10, 10), sq(150, 10, 10)])
        );
        z.verify().unwrap();
    }

    #[test]
    fn tags_survive_incremental_updates() {
        let zone = TaggedRing::uniform(sq(0, 0, 100), 7);
        let mut z = ZoneFill::new(&zone, FillRule::NonZero).unwrap();
        z.set_auto_rebuild(false);
        let a = TaggedRing::uniform(sq(-5, 20, 30), 11);
        let b = TaggedRing::uniform(sq(50, 50, 20), 12);
        z.insert(1, &a).unwrap();
        z.insert(2, &b).unwrap();
        let want = Boolean::new()
            .subject(&zone, FillRule::NonZero)
            .clip(&[a.clone(), b], FillRule::NonZero)
            .op(Op::Difference)
            .execute_tagged()
            .unwrap();
        assert_eq!(z.fill_tagged(), want);
        z.remove(2);
        let want = Boolean::new()
            .subject(&zone, FillRule::NonZero)
            .clip(&a, FillRule::NonZero)
            .op(Op::Difference)
            .execute_tagged()
            .unwrap();
        assert_eq!(z.fill_tagged(), want);
        z.verify().unwrap();
    }

    #[test]
    fn touching_obstacles_and_pinches() {
        // Obstacles meeting at corners create pinch vertices in the fill.
        let zone = sq(0, 0, 40);
        let mut z = ZoneFill::new(&zone, FillRule::NonZero).unwrap();
        z.set_auto_rebuild(false);
        let shapes = [
            sq(10, 10, 10),
            sq(20, 20, 10),
            sq(0, 20, 10),
            sq(20, 0, 10),
            sq(10, 30, 10),
            sq(30, 10, 10),
        ];
        for (i, s) in shapes.iter().enumerate() {
            z.insert(i as u64, s).unwrap();
            assert_eq!(z.fill_tree(), reference_tree(&zone, &shapes[..=i]));
            z.verify().unwrap();
        }
        for i in (0..shapes.len()).rev() {
            z.remove(i as u64);
            assert_eq!(z.fill_tree(), reference_tree(&zone, &shapes[..i]));
            z.verify().unwrap();
        }
        assert_eq!(z.rebuild_count(), 1);
    }

    #[test]
    fn rounding_cascade_through_a_chain() {
        // A fan of nearly parallel slivers crossed by a cutter: rounding at the crossings
        // reroutes neighbours, which activates more pixels, and so on.
        let zone = sq(0, 0, 1000);
        let mut z = ZoneFill::new(&zone, FillRule::NonZero).unwrap();
        z.set_auto_rebuild(false);
        let mut obst = Vec::new();
        for k in 0..12i64 {
            let r = Ring::from([
                (100 + k, 100),
                (900 + k, 701 + 2 * k),
                (900 + k, 703 + 2 * k),
            ]);
            z.insert(k as u64, &r).unwrap();
            obst.push(r);
        }
        assert_eq!(z.fill(), reference(&zone, &obst));
        let cut = Ring::from([(500, 0), (520, 1000), (521, 1000)]);
        z.insert(100, &cut).unwrap();
        let mut with = obst.clone();
        with.push(cut);
        assert_eq!(z.fill_tree(), reference_tree(&zone, &with));
        z.verify().unwrap();
        z.remove(100);
        assert_eq!(z.fill_tree(), reference_tree(&zone, &obst));
        z.verify().unwrap();
        assert_eq!(z.rebuild_count(), 1);
    }

    #[test]
    fn far_and_huge_obstacles() {
        // Long segments (grid long list), cells outside the dense extent, long rays.
        let zone = sq(0, 0, 50);
        let mut z = ZoneFill::new(&zone, FillRule::NonZero).unwrap();
        z.set_auto_rebuild(false);
        let huge = Ring::from([(-(1 << 39), 20), (1 << 39, 21), (1 << 39, 25)]);
        let far = sq(1 << 38, 1 << 38, 1000);
        let small = sq(10, 10, 5);
        z.insert(1, &huge).unwrap();
        z.insert(2, &far).unwrap();
        z.insert(3, &small).unwrap();
        let all = [huge.clone(), far.clone(), small.clone()];
        assert_eq!(z.fill_tree(), reference_tree(&zone, &all));
        z.verify().unwrap();
        z.remove(1);
        assert_eq!(z.fill_tree(), reference_tree(&zone, &[far, small]));
        z.verify().unwrap();
    }

    #[test]
    fn grid_cells_are_conservative() {
        // Every cell whose grown square meets the segment is visited.
        let p = Point::new;
        for &(a, b) in &[
            (p(0, 0), p(1000, 377)),
            (p(-50, 900), p(800, -3)),
            (p(5, 5), p(5, 900)),
            (p(-7, 3), p(2000, 3)),
            (p(3, 1), p(2, 999)),
        ] {
            for sh in [0u32, 3, 6, 9] {
                let mut cells = Vec::new();
                for_seg_cells(a, b, sh, &mut |c| cells.push(c));
                let (x0, x1) = ((a.x.min(b.x) - 2) >> sh, (a.x.max(b.x) + 2) >> sh);
                let (y0, y1) = ((a.y.min(b.y) - 2) >> sh, (a.y.max(b.y) + 2) >> sh);
                for cy in y0..=y1 {
                    for cx in x0..=x1 {
                        if segment_meets_rect(a, b, &cell_rect_grown((cx, cy), sh)) {
                            assert!(
                                cells.contains(&(cx, cy)),
                                "{a:?}-{b:?} sh {sh} misses {cx},{cy}"
                            );
                        }
                    }
                }
            }
        }
    }
}
