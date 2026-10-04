//! Prepared geometries: a spatial index over the segments and rings of a geometry that is
//! queried many times (a zone fill checked against every pad and track in DRC).
//!
//! Every query gives exactly the result of the corresponding free function, closest points
//! included; only the work changes. The free functions read every segment of both operands
//! on every call (to test ranges, collect the segments near the other operand, and locate
//! points by walking every ring). A [`Prepared`] geometry does that once:
//!
//! * the segments near the other operand come from a bounding-volume hierarchy, in the
//!   same order as a full scan would list them, and are handed to the same pair search;
//! * points are located with a ray towards the nearest side of the bounding box, which only
//!   meets the segments it crosses; winding numbers and even-odd parities are properties of
//!   the point and the closed rings, so any ray gives the counts the full walk gives;
//! * [`contains`](Prepared::contains) builds the arrangement of the segments near the other
//!   operand only, and takes the parity of every face that matters from a ray as well;
//! * [`distance`](Prepared::distance) finds the minimum by branch-and-bound over the index,
//!   then reproduces which pair of closest points the full search reports (the first one in
//!   its scan order).

use crate::arrangement::{Arrangement, InEdge, Noding};
use crate::dir::Dir;
use crate::distance::{Bvh, Closest, SqDist, rect_gap2, segment_segment, segment_segment_lt};
use crate::geom::{Path, Point, PointF, PolyTree, Polygon, Rect, Ring};
use crate::predicates::{on_segment, orient, segments_intersect};
use crate::query::{
    EMPTY_RECT, Geometry, Location, Segment, any_pair, collect_segments, dir_sweep, in_range,
};
use std::sync::OnceLock;

pub(crate) mod sealed {
    /// How the region of a geometry is defined by its rings.
    pub enum Mode {
        /// Points and polylines: no region.
        Linear,
        /// Polygons (outer ring, then holes) in `rings`, grouped by `polys`.
        Polygons,
        /// Every ring counts (odd number of enclosing rings means inside).
        Tree,
    }

    /// Ring structure of a geometry, in [`visit_segments`](crate::Geometry::visit_segments)
    /// order.
    pub struct Layout {
        pub mode: Mode,
        /// Vertex count of every ring.
        pub rings: Vec<u32>,
        /// Ring count of every polygon ([`Mode::Polygons`]).
        pub polys: Vec<u32>,
    }

    pub trait Sealed {
        fn layout(&self) -> Layout;
    }
}

use sealed::{Layout, Mode};

/// Geometries that can be [`Prepared`]: [`Ring`], [`Polygon`], polygon sets, [`PolyTree`],
/// [`Path`], [`Segment`] and [`Point`]. Sealed.
pub trait Preparable: Geometry + sealed::Sealed {}

impl sealed::Sealed for Ring {
    fn layout(&self) -> Layout {
        Layout {
            mode: Mode::Polygons,
            rings: vec![self.0.len() as u32],
            polys: vec![1],
        }
    }
}

impl sealed::Sealed for Polygon {
    fn layout(&self) -> Layout {
        Layout {
            mode: Mode::Polygons,
            rings: self.rings().map(|r| r.0.len() as u32).collect(),
            polys: vec![1 + self.holes.len() as u32],
        }
    }
}

impl sealed::Sealed for [Polygon] {
    fn layout(&self) -> Layout {
        Layout {
            mode: Mode::Polygons,
            rings: self
                .iter()
                .flat_map(|p| p.rings().map(|r| r.0.len() as u32))
                .collect(),
            polys: self.iter().map(|p| 1 + p.holes.len() as u32).collect(),
        }
    }
}

impl sealed::Sealed for Vec<Polygon> {
    fn layout(&self) -> Layout {
        self.as_slice().layout()
    }
}

impl sealed::Sealed for PolyTree {
    fn layout(&self) -> Layout {
        Layout {
            mode: Mode::Tree,
            rings: self.nodes.iter().map(|n| n.ring.0.len() as u32).collect(),
            polys: Vec::new(),
        }
    }
}

macro_rules! linear_layout {
    ($($t:ty),*) => {$(
        impl sealed::Sealed for $t {
            fn layout(&self) -> Layout {
                Layout { mode: Mode::Linear, rings: Vec::new(), polys: Vec::new() }
            }
        }
        impl Preparable for $t {}
    )*};
}

linear_layout!(Path, Segment, Point);

impl<T: sealed::Sealed + ?Sized> sealed::Sealed for &T {
    fn layout(&self) -> Layout {
        (**self).layout()
    }
}

impl Preparable for Ring {}
impl Preparable for Polygon {}
impl Preparable for [Polygon] {}
impl Preparable for Vec<Polygon> {}
impl Preparable for PolyTree {}
impl<T: Preparable + ?Sized> Preparable for &T {}

/// A segment with its bounding box.
type Seg = (Point, Point, Rect);

/// Below this many segments (both operands together) the free functions are called
/// directly: they are fast there, and [`distance`](crate::distance) orders its search
/// differently for small inputs.
const SMALL: usize = 512;

/// Bounding-volume hierarchy over item indices: interior nodes have `count == 0` and
/// children at `first`, `first + 1`; leaves cover `items[first..first + count]`.
pub(crate) struct Index {
    pub nodes: Vec<(Rect, u32, u32)>,
    pub items: Vec<u32>,
}

const LEAF: usize = 8;

impl Index {
    fn build(segs: &[Seg]) -> Index {
        Self::of_boxes(segs.len(), |i| segs[i].2)
    }

    /// Builds the hierarchy over items `0..n` with boxes `bx(i)`.
    pub(crate) fn of_boxes(n: usize, bx: impl Fn(usize) -> Rect) -> Index {
        let mut it: Vec<(Rect, u32)> = (0..n).map(|i| (bx(i), i as u32)).collect();
        let union = |s: &[(Rect, u32)]| s.iter().fold(EMPTY_RECT, |a, b| a.union(&b.0));
        let mut nodes = vec![(union(&it), 0u32, n as u32)];
        let mut stack = vec![0usize];
        while let Some(ni) = stack.pop() {
            let (b, first, count) = nodes[ni];
            let (first, count) = (first as usize, count as usize);
            if count <= LEAF {
                continue;
            }
            let part = &mut it[first..first + count];
            let mid = count / 2;
            // (Wide arithmetic: coordinates may be out of range here.)
            if b.max.x as i128 - b.min.x as i128 >= b.max.y as i128 - b.min.y as i128 {
                part.select_nth_unstable_by_key(mid, |e| e.0.min.x as i128 + e.0.max.x as i128);
            } else {
                part.select_nth_unstable_by_key(mid, |e| e.0.min.y as i128 + e.0.max.y as i128);
            }
            let (l, rr) = (union(&part[..mid]), union(&part[mid..]));
            let c = nodes.len() as u32;
            nodes.push((l, first as u32, mid as u32));
            nodes.push((rr, (first + mid) as u32, (count - mid) as u32));
            nodes[ni] = (b, c, 0);
            stack.push(c as usize);
            stack.push(c as usize + 1);
        }
        Index {
            nodes,
            items: it.into_iter().map(|e| e.1).collect(),
        }
    }

    /// Calls `f` for every item in a leaf whose box passes `node` (tested on every node on
    /// the way down).
    pub(crate) fn visit(&self, node: impl Fn(&Rect) -> bool, mut f: impl FnMut(u32)) {
        if self.items.is_empty() {
            return;
        }
        let mut stack = vec![0u32];
        while let Some(ni) = stack.pop() {
            let (b, first, count) = self.nodes[ni as usize];
            if !node(&b) {
                continue;
            }
            if count == 0 {
                stack.push(first + 1);
                stack.push(first);
            } else {
                for &i in &self.items[first as usize..(first + count) as usize] {
                    f(i);
                }
            }
        }
    }
}

/// Projected extent of a segment set along a direction: (lowest start, highest end, total
/// width), as [`crate::dir::separating`] computes it.
type Ext = (i128, i128, i128);

fn ext_of(d: Dir, segs: impl Iterator<Item = (Point, Point)>) -> Ext {
    let (mut lo, mut hi, mut w) = (i128::MAX, i128::MIN, 0i128);
    for s in segs {
        let (l, h) = d.range(&s);
        lo = lo.min(l);
        hi = hi.max(h);
        w += h - l;
    }
    (lo, hi, w)
}

/// The normal of a segment, reduced (as [`Dir::best_of`] picks it).
fn normal(s: (Point, Point)) -> Dir {
    let (dx, dy) = (s.1.x - s.0.x, s.1.y - s.0.y);
    let (mut a, mut b) = (dx.unsigned_abs(), dy.unsigned_abs());
    while b != 0 {
        (a, b) = (b, a % b);
    }
    let g = a.max(1) as i64;
    Dir {
        nx: -dy / g,
        ny: dx / g,
    }
}

/// Aggregates of all segments for reproducing [`crate::dir::separating`] without a scan.
struct DirAgg {
    fixed: [Ext; 4],
    /// The last longest segment's squared length, normal and extent along it.
    longest: (i128, Dir, Ext),
}

/// Rotations taking a ray direction to +x (right, left, up, down), as
/// `(x, y) -> (c * x - s * y, s * x + c * y)`.
const ROT: [(i64, i64); 4] = [(1, 0), (-1, 0), (0, -1), (0, 1)];

#[inline]
fn rot(p: Point, r: (i64, i64)) -> Point {
    let (c, s) = r;
    Point::new(c * p.x - s * p.y, s * p.x + c * p.y)
}

/// A geometry prepared for repeated queries: a spatial index over its segments and rings.
///
/// Build it once for a geometry queried many times (a zone fill checked against every pad,
/// track and via in DRC); every method returns exactly what the free function of the same
/// name returns with the prepared geometry as first argument, closest points included, at
/// a cost that follows the size of the other operand and of the prepared geometry near it
/// rather than the size of the whole prepared geometry.
///
/// Building costs about as much as a few free-function queries (`O(n log n)` for `n`
/// segments); [`contains`](Self::contains) and [`distance`](Self::distance) also compute a
/// little more on their first call. A `Prepared` is `Send + Sync` when the geometry is.
///
/// ```
/// use polyclip::*;
/// let zone = Polygon::new(
///     Ring::from([(0, 0), (1000, 0), (1000, 1000), (0, 1000)]),
///     vec![Ring::from([(400, 400), (400, 600), (600, 600), (600, 400)])],
/// );
/// let prepared = Prepared::new(&zone);
/// let pad = Ring::from([(450, 450), (550, 450), (550, 550), (450, 550)]);
/// assert!(!prepared.intersects(&pad));
/// assert!(prepared.distance_less_than(&pad, 51));
/// assert!(!prepared.distance_less_than(&pad, 50));
/// assert_eq!(prepared.locate(Point::new(500, 500)), Location::Outside);
/// assert_eq!(prepared.distance_less_than(&pad, 51), distance_less_than(&zone, &pad, 51));
/// ```
pub struct Prepared<'a, G: Preparable + ?Sized> {
    geom: &'a G,
    bbox: Option<Rect>,
    /// Union of all segment boxes (holes and empty outer rings included).
    full: Option<Rect>,
    in_range: bool,
    areal: bool,
    segs: Vec<Seg>,
    index: Index,
    /// Ring of every segment (when the layout matched the segments).
    ring_of: Vec<u32>,
    mode: Mode,
    /// Polygon of every ring, and first ring of every polygon ([`Mode::Polygons`]).
    poly_of: Vec<u32>,
    poly_first: Vec<u32>,
    /// One point per connected component, in [`Geometry::component_points`] order.
    comps: Vec<Point>,
    /// Whether two segments cross properly (computed on demand).
    crossing: OnceLock<bool>,
    dir_agg: OnceLock<DirAgg>,
}

impl<'a, G: Preparable + ?Sized> Prepared<'a, G> {
    /// Indexes `geom`.
    pub fn new(geom: &'a G) -> Self {
        let mut segs: Vec<Seg> = Vec::new();
        geom.visit_segments(&mut |a, b| segs.push((a, b, Rect::new(a, b))));
        let in_range = segs.iter().all(|s| s.0.in_range() && s.1.in_range());
        let full = segs.iter().map(|s| s.2).reduce(|a, b| a.union(&b));
        let layout = geom.layout();
        let mut mode = layout.mode;
        let mut ring_of = Vec::new();
        let (mut poly_of, mut poly_first) = (Vec::new(), Vec::new());
        if !matches!(mode, Mode::Linear) {
            ring_of.reserve(segs.len());
            for (r, &n) in layout.rings.iter().enumerate() {
                ring_of.extend(core::iter::repeat_n(r as u32, n as usize));
            }
            let rings_ok = ring_of.len() == segs.len();
            let polys_ok = !matches!(mode, Mode::Polygons)
                || layout.polys.iter().map(|&n| n as usize).sum::<usize>() == layout.rings.len();
            if rings_ok && polys_ok {
                for (p, &n) in layout.polys.iter().enumerate() {
                    poly_first.push(poly_of.len() as u32);
                    poly_of.extend(core::iter::repeat_n(p as u32, n as usize));
                }
            } else {
                // Not expected for the sealed implementations; locate through the geometry.
                ring_of.clear();
                mode = Mode::Linear;
            }
        }
        let mut comps = Vec::new();
        geom.component_points(&mut |p| comps.push(p));
        let index = Index::build(&segs);
        Prepared {
            geom,
            bbox: geom.bbox(),
            full,
            in_range,
            areal: geom.is_areal(),
            segs,
            index,
            ring_of,
            mode,
            poly_of,
            poly_first,
            comps,
            crossing: OnceLock::new(),
            dir_agg: OnceLock::new(),
        }
    }

    /// The prepared geometry.
    pub fn geometry(&self) -> &'a G {
        self.geom
    }

    /// Segments whose box meets `r`, in visiting order (as [`collect_segments`] lists them).
    fn collect(&self, r: &Rect) -> Vec<Seg> {
        let mut idx: Vec<u32> = Vec::new();
        self.index.visit(
            |b| b.intersects(r),
            |i| {
                if self.segs[i as usize].2.intersects(r) {
                    idx.push(i)
                }
            },
        );
        idx.sort_unstable();
        idx.iter().map(|&i| self.segs[i as usize]).collect()
    }

    /// The ray direction (index into [`ROT`]) leaving the segments' box soonest from `p`.
    fn ray(&self, p: Point) -> usize {
        let Some(f) = self.full else { return 0 };
        let d = [
            f.max.x.saturating_sub(p.x),
            p.x.saturating_sub(f.min.x),
            f.max.y.saturating_sub(p.y),
            p.y.saturating_sub(f.min.y),
        ];
        (0..4).min_by_key(|&k| d[k]).unwrap_or(0)
    }

    /// Calls `f(segment)` for every segment whose box meets the ray from `p` in direction
    /// `k` (`p` included).
    fn ray_segments(&self, p: Point, k: usize, mut f: impl FnMut(u32)) {
        let Some(full) = self.full else { return };
        let r = match k {
            0 => Rect::new(p, Point::new(full.max.x.max(p.x), p.y)),
            1 => Rect::new(Point::new(full.min.x.min(p.x), p.y), p),
            2 => Rect::new(p, Point::new(p.x, full.max.y.max(p.y))),
            _ => Rect::new(Point::new(p.x, full.min.y.min(p.y)), p),
        };
        self.index.visit(
            |b| b.intersects(&r),
            |i| {
                if self.segs[i as usize].2.intersects(&r) {
                    f(i)
                }
            },
        );
    }

    /// `g.locate(p)` for the prepared geometry `g` (coordinates in range).
    fn locate_raw(&self, p: Point) -> Location {
        if matches!(self.mode, Mode::Linear) {
            if self.areal {
                // Only if the layout did not match the segments.
                return self.geom.locate(p);
            }
            let mut on = false;
            let r = Rect::new(p, p);
            self.index.visit(
                |b| b.intersects(&r),
                |i| {
                    let s = &self.segs[i as usize];
                    on |= s.2.contains_point(p) && on_segment(s.0, s.1, p);
                },
            );
            return if on {
                Location::OnBoundary
            } else {
                Location::Outside
            };
        }
        // Winding number of every ring the ray meets, with `ring_winding`'s rule applied in
        // a frame where the ray points along +x.
        let k = self.ray(p);
        let rp = rot(p, ROT[k]);
        let mut hits: Vec<(u32, i32, bool)> = Vec::new();
        self.ray_segments(p, k, |i| {
            let s = &self.segs[i as usize];
            let (a, b) = (rot(s.0, ROT[k]), rot(s.1, ROT[k]));
            let (mut w, mut on) = (0, false);
            if a.y <= rp.y {
                if b.y > rp.y {
                    let o = orient(a, b, rp);
                    if o > 0 {
                        w = 1;
                    } else if o == 0 {
                        on = true;
                    }
                } else if b.y == rp.y && on_segment(a, b, rp) {
                    on = true;
                }
            } else if b.y <= rp.y {
                let o = orient(a, b, rp);
                if o < 0 {
                    w = -1;
                } else if o == 0 {
                    on = true;
                }
            }
            if w != 0 || on {
                hits.push((self.ring_of[i as usize], w, on));
            }
        });
        hits.sort_unstable_by_key(|h| h.0);
        // Per ring: (ring, winding, on).
        let mut rings: Vec<(u32, i32, bool)> = Vec::new();
        for (r, w, on) in hits {
            match rings.last_mut() {
                Some(l) if l.0 == r => {
                    l.1 += w;
                    l.2 |= on;
                }
                _ => rings.push((r, w, on)),
            }
        }
        let status = |&(_, w, on): &(u32, i32, bool)| {
            if on {
                Location::OnBoundary
            } else if w != 0 {
                Location::Inside
            } else {
                Location::Outside
            }
        };
        match self.mode {
            Mode::Tree => {
                if rings.iter().any(|r| r.2) {
                    Location::OnBoundary
                } else if rings.iter().filter(|r| r.1 != 0).count() % 2 == 1 {
                    Location::Inside
                } else {
                    Location::Outside
                }
            }
            _ => {
                // `locate` of a polygon list: inside one polygon, else on one's boundary.
                let mut res = Location::Outside;
                let mut i = 0;
                while i < rings.len() {
                    let poly = self.poly_of[rings[i].0 as usize];
                    let mut j = i + 1;
                    while j < rings.len() && self.poly_of[rings[j].0 as usize] == poly {
                        j += 1;
                    }
                    let group = &rings[i..j];
                    i = j;
                    // The outer ring decides first, then the first hole that is not outside.
                    if group[0].0 != self.poly_first[poly as usize] {
                        continue;
                    }
                    let here = match status(&group[0]) {
                        Location::Inside => group[1..]
                            .iter()
                            .map(status)
                            .find(|s| *s != Location::Outside)
                            .map_or(Location::Inside, |s| match s {
                                Location::Inside => Location::Outside,
                                other => other,
                            }),
                        other => other,
                    };
                    match here {
                        Location::Inside => return Location::Inside,
                        Location::OnBoundary => res = Location::OnBoundary,
                        Location::Outside => {}
                    }
                }
                res
            }
        }
    }

    /// Location of `p` relative to the prepared geometry: [`locate`](crate::locate)`(g, p)`.
    pub fn locate(&self, p: Point) -> Location {
        if !p.in_range() || !self.in_range {
            return Location::Outside;
        }
        self.locate_raw(p)
    }

    /// `any_component_inside(b, g)`: a component point of `g` in the region `b`.
    fn component_in<B: Geometry + ?Sized>(&self, b: &B) -> Option<Point> {
        if !b.is_areal() {
            return None;
        }
        let ob = b.bbox()?;
        let pts: Vec<Point> = self
            .comps
            .iter()
            .copied()
            .filter(|&p| ob.contains_point(p))
            .collect();
        b.first_not_outside(&pts)
    }

    /// `any_component_inside(g, b)`: a component point of `b` in the prepared region.
    fn first_inside<B: Geometry + ?Sized>(&self, b: &B) -> Option<Point> {
        if !self.areal {
            return None;
        }
        let ob = self.bbox?;
        let mut pts = Vec::new();
        b.component_points(&mut |p| {
            if ob.contains_point(p) {
                pts.push(p);
            }
        });
        pts.into_iter()
            .find(|&p| self.locate_raw(p) != Location::Outside)
    }

    /// `true` when the prepared geometry and `b` share a point:
    /// [`intersects`](crate::intersects)`(g, b)`.
    pub fn intersects<B: Geometry + ?Sized>(&self, b: &B) -> bool {
        let (Some(ba), Some(bb)) = (self.bbox, b.bbox()) else {
            return false;
        };
        if ![ba.min, ba.max, bb.min, bb.max]
            .iter()
            .all(|p| p.in_range())
        {
            return false;
        }
        if !ba.intersects(&bb) || !self.in_range || !in_range(b) {
            return false;
        }
        let mut sa = self.collect(&bb);
        let mut sb = collect_segments(b, &ba);
        if any_pair(&mut sa, &mut sb, 0, |x, y| {
            segments_intersect(x.0, x.1, y.0, y.1)
        }) {
            return true;
        }
        self.component_in(b).is_some() || self.first_inside(b).is_some()
    }

    /// `true` when the distance between the prepared geometry and `b` is less than `d`:
    /// [`distance_less_than`](crate::distance_less_than)`(g, b, d)`.
    pub fn distance_less_than<B: Geometry + ?Sized>(&self, b: &B, d: i64) -> bool {
        if d <= 0 {
            return false;
        }
        let (Some(ba), Some(bb)) = (self.bbox, b.bbox()) else {
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
        if !self.in_range || !in_range(b) {
            return false;
        }
        let mut sa = self.collect(&bb.expand(d));
        let mut sb = collect_segments(b, &ba.expand(d));
        if any_pair(&mut sa, &mut sb, d, |x, y| {
            rect_gap2(&x.2, &y.2) < d2 && segment_segment_lt(x.0, x.1, y.0, y.1, d2)
        }) {
            return true;
        }
        self.component_in(b).is_some() || self.first_inside(b).is_some()
    }

    /// Whether two segments of the prepared geometry cross properly (then no arrangement
    /// with exact noding exists, and [`contains`](crate::contains) is `false`).
    fn self_crossing(&self) -> bool {
        *self.crossing.get_or_init(|| {
            let segs: Vec<(Point, Point)> = self
                .segs
                .iter()
                .filter(|s| s.0 != s.1)
                .map(|s| (s.0, s.1))
                .collect();
            crate::node::node_exact(&segs).is_err()
        })
    }

    /// Even-odd parity of the prepared rings just below the arrangement edge `lo-hi`
    /// (`lo < hi`; right of it when vertical, as the sweep sees it), counted along a
    /// vertical (horizontal) ray from the middle of its first unit of length.
    fn parity_below(&self, lo: Point, hi: Point) -> bool {
        let Some(full) = self.full else { return false };
        let mut odd = false;
        if lo.x < hi.x {
            // At x* = lo.x + 1/2: count segments strictly below the edge there (downwards)
            // or on or above it (upwards). No segment is vertical at x*.
            let (dxe, dye) = ((hi.x - lo.x) as i128, (hi.y - lo.y) as i128);
            let up = full.max.y - lo.y.max(hi.y) < lo.y.min(hi.y) - full.min.y;
            let (ylo, yhi) = (lo.y.min(hi.y), lo.y.max(hi.y));
            self.index.visit(
                |b| {
                    b.min.x <= lo.x
                        && b.max.x > lo.x
                        && if up { b.max.y >= ylo } else { b.min.y <= yhi }
                },
                |i| {
                    let s = &self.segs[i as usize];
                    let (u, v) = if s.0.x < s.1.x {
                        (s.0, s.1)
                    } else {
                        (s.1, s.0)
                    };
                    if !(u.x <= lo.x && v.x > lo.x) {
                        return;
                    }
                    let (dxf, dyf) = ((v.x - u.x) as i128, (v.y - u.y) as i128);
                    let t = 2 * (lo.x - u.x) as i128 + 1;
                    // sign of y_f(x*) - y_e(x*), times 2 * dxf * dxe > 0.
                    let c = 2 * dxf * dxe * (u.y - lo.y) as i128 + dyf * t * dxe - dye * dxf;
                    if (c < 0) != up {
                        odd = !odd;
                    }
                },
            );
        } else {
            // Vertical edge at x = c: right side, at y* = lo.y + 1/2: count segments
            // strictly right of it (rightwards) or on or left of it (leftwards).
            let c = lo.x;
            let left = c - full.min.x < full.max.x - c;
            self.index.visit(
                |b| {
                    b.min.y <= lo.y
                        && b.max.y > lo.y
                        && if left { b.min.x <= c } else { b.max.x >= c }
                },
                |i| {
                    let s = &self.segs[i as usize];
                    let (u, v) = if s.0.y < s.1.y {
                        (s.0, s.1)
                    } else {
                        (s.1, s.0)
                    };
                    if !(u.y <= lo.y && v.y > lo.y) {
                        return;
                    }
                    let (dxf, dyf) = ((v.x - u.x) as i128, (v.y - u.y) as i128);
                    let t = 2 * (lo.y - u.y) as i128 + 1;
                    // sign of x_f(y*) - c, times 2 * dyf > 0.
                    let k = 2 * dyf * (u.x - c) as i128 + dxf * t;
                    if (k > 0) != left {
                        odd = !odd;
                    }
                },
            );
        }
        odd
    }

    /// Even-odd location of `m2 / 2` (doubled coordinates) relative to all segments:
    /// `locate_doubled` over every segment of the prepared geometry.
    fn locate_doubled(&self, m2: Point) -> Location {
        let Some(full) = self.full else {
            return Location::Outside;
        };
        let f2 = Rect {
            min: Point::new(2 * full.min.x, 2 * full.min.y),
            max: Point::new(2 * full.max.x, 2 * full.max.y),
        };
        let d = [
            f2.max.x.saturating_sub(m2.x),
            m2.x.saturating_sub(f2.min.x),
            f2.max.y.saturating_sub(m2.y),
            m2.y.saturating_sub(f2.min.y),
        ];
        let k = (0..4).min_by_key(|&k| d[k]).unwrap_or(0);
        // The ray from m2 / 2, in undoubled coordinates rounded outwards.
        let (lo, hi) = (
            Point::new(m2.x.div_euclid(2), m2.y.div_euclid(2)),
            Point::new((m2.x + 1).div_euclid(2), (m2.y + 1).div_euclid(2)),
        );
        let r = match k {
            0 => Rect::new(lo, Point::new(full.max.x.max(hi.x), hi.y)),
            1 => Rect::new(Point::new(full.min.x.min(lo.x), lo.y), hi),
            2 => Rect::new(lo, Point::new(hi.x, full.max.y.max(hi.y))),
            _ => Rect::new(Point::new(lo.x, full.min.y.min(lo.y)), hi),
        };
        let rm = rot(m2, ROT[k]);
        let (mut inside, mut on) = (false, false);
        self.index.visit(
            |b| b.intersects(&r),
            |i| {
                let s = &self.segs[i as usize];
                if !s.2.intersects(&r) {
                    return;
                }
                let a = rot(Point::new(2 * s.0.x, 2 * s.0.y), ROT[k]);
                let b = rot(Point::new(2 * s.1.x, 2 * s.1.y), ROT[k]);
                if on_segment(a, b, rm) {
                    on = true;
                } else if (a.y > rm.y) != (b.y > rm.y) {
                    let o = orient(a, b, rm);
                    if (o > 0) == (b.y > a.y) {
                        inside = !inside;
                    }
                }
            },
        );
        if on {
            Location::OnBoundary
        } else if inside {
            Location::Inside
        } else {
            Location::Outside
        }
    }

    /// `true` when every point of `b` belongs to the prepared geometry:
    /// [`contains`](crate::contains)`(g, b)`.
    pub fn contains<B: Geometry + ?Sized>(&self, b: &B) -> bool {
        let Some(bb) = b.bbox() else {
            return true;
        };
        let Some(ba) = self.bbox else {
            return false;
        };
        if ![ba.min, ba.max, bb.min, bb.max]
            .iter()
            .all(|p| p.in_range())
            || !ba.contains_rect(&bb)
            || !self.in_range
            || !in_range(b)
        {
            return false;
        }
        if !self.areal {
            let mut pts = Vec::new();
            let mut linear = false;
            b.visit_segments(&mut |p, q| {
                if p == q { pts.push(p) } else { linear = true }
            });
            return !linear
                && !b.is_areal()
                && pts.iter().all(|&p| self.locate_raw(p) != Location::Outside);
        }
        if !b.is_areal() {
            return self.linear_inside(b);
        }
        let mut isolated = Vec::new();
        let mut edges: Vec<InEdge> = Vec::new();
        // All of `b`'s segments (holes may stick out of the outer ring's box).
        let mut sbox = bb;
        b.visit_segments(&mut |p, q| {
            sbox.add_point(p);
            sbox.add_point(q);
        });
        for s in self.collect(&sbox) {
            if s.0 != s.1 {
                edges.push(InEdge {
                    a: s.0,
                    b: s.1,
                    tag: 0,
                    operand: 0,
                });
            }
        }
        b.visit_segments(&mut |p, q| {
            if p == q {
                isolated.push(p);
            } else {
                edges.push(InEdge {
                    a: p,
                    b: q,
                    tag: 0,
                    operand: 1,
                });
            }
        });
        if isolated
            .iter()
            .any(|&p| self.locate_raw(p) == Location::Outside)
        {
            return false;
        }
        // The full arrangement fails on any proper crossing: among the segments near `b`
        // (found here) or anywhere in the prepared geometry.
        let Ok(arr) = Arrangement::build(&edges, Noding::Exact) else {
            return false;
        };
        if self.self_crossing() {
            return false;
        }
        // Edges with `b` on a side lie in `b`'s box, where this arrangement and the full one
        // agree (every vertex or edge meeting them is near `b`), and so do `b`'s windings
        // (all of `b` is here). The parity of the prepared rings comes from a ray instead.
        for k in 0..arr.edges.len() {
            let (wb, wa) = arr.sides(k);
            let (ob, oa) = (wb[1] & 1 != 0, wa[1] & 1 != 0);
            if !ob && !oa {
                continue;
            }
            let e = &arr.edges[k];
            let pb = self.parity_below(e.lo, e.hi);
            let pa = pb ^ (e.delta[0] & 1 != 0);
            if (ob && !pb) || (oa && !pa) {
                return false;
            }
        }
        true
    }

    /// `linear_inside(g, b)` with the segments near `b` only.
    fn linear_inside<B: Geometry + ?Sized>(&self, b: &B) -> bool {
        let mut sb: Vec<Seg> = Vec::new();
        let mut isolated: Vec<Point> = Vec::new();
        let mut sbox = EMPTY_RECT;
        b.visit_segments(&mut |p, q| {
            sbox.add_point(p);
            sbox.add_point(q);
            if p == q {
                isolated.push(p);
            } else {
                sb.push((p, q, Rect::new(p, q)));
            }
        });
        if isolated
            .iter()
            .any(|&p| self.locate_raw(p) == Location::Outside)
        {
            return false;
        }
        // Pairs need meeting boxes: only segments near `b` take part.
        let mut sa = self.collect(&sbox);
        let mut touches: Vec<((Point, Point), i128, Point)> = Vec::new();
        let crossed = any_pair(&mut sa, &mut sb.clone(), 0, |x, y| {
            if crate::predicates::segments_cross_properly(x.0, x.1, y.0, y.1) {
                return true;
            }
            for v in [x.0, x.1] {
                if crate::predicates::in_segment_interior(y.0, y.1, v) {
                    touches.push(((y.0, y.1), crate::predicates::dist2(y.0, v), v));
                }
            }
            false
        });
        if crossed {
            return false;
        }
        touches.sort_unstable();
        touches.dedup();
        let mut t = 0usize;
        sb.sort_unstable_by_key(|s| (s.0, s.1));
        for &(p, q, _) in &sb {
            while t < touches.len() && touches[t].0 < (p, q) {
                t += 1;
            }
            let mut cur = p;
            let mut pieces: Vec<(Point, Point)> = Vec::new();
            while t < touches.len() && touches[t].0 == (p, q) {
                let v = touches[t].2;
                if v != cur {
                    pieces.push((cur, v));
                    cur = v;
                }
                t += 1;
            }
            pieces.push((cur, q));
            for (u, w) in pieces {
                let m2 = Point::new(u.x + w.x, u.y + w.y);
                if self.locate_doubled(m2) == Location::Outside {
                    return false;
                }
            }
        }
        true
    }

    fn dir_agg(&self) -> &DirAgg {
        self.dir_agg.get_or_init(|| {
            let it = || self.segs.iter().map(|s| (s.0, s.1));
            let fixed = Dir::FIXED.map(|d| ext_of(d, it()));
            // `max_by_key` keeps the last maximum.
            let mut best: Option<(i128, usize)> = None;
            for (i, s) in self.segs.iter().enumerate() {
                let l = crate::predicates::dist2(s.0, s.1);
                if best.is_none_or(|(b, _)| l >= b) {
                    best = Some((l, i));
                }
            }
            let (l, i) = best.unwrap_or((0, 0));
            let n = self
                .segs
                .get(i)
                .map_or(Dir { nx: 0, ny: 0 }, |s| normal((s.0, s.1)));
            DirAgg {
                fixed,
                longest: (l, n, ext_of(n, it())),
            }
        })
    }

    /// `crate::dir::separating(all segments, sb)` from the aggregates, or `None` when it
    /// needs the extent along a direction not aggregated.
    fn separating(&self, sb: &[Seg]) -> Option<Dir> {
        let agg = self.dir_agg();
        let known = |d: Dir| -> Option<Ext> {
            let neg = |e: Ext| (-e.1, -e.0, e.2);
            for (k, f) in Dir::FIXED.iter().enumerate() {
                if *f == d {
                    return Some(agg.fixed[k]);
                }
                if f.nx == -d.nx && f.ny == -d.ny {
                    return Some(neg(agg.fixed[k]));
                }
            }
            let n = agg.longest.1;
            if n == d {
                Some(agg.longest.2)
            } else if n.nx == -d.nx && n.ny == -d.ny {
                Some(neg(agg.longest.2))
            } else {
                None
            }
        };
        let itb = || sb.iter().map(|s| (s.0, s.1));
        // Dir::best_of over both: the last longest segment's normal as extra candidate.
        let mut lb: Option<(i128, usize)> = None;
        for (i, s) in sb.iter().enumerate() {
            let l = crate::predicates::dist2(s.0, s.1);
            if lb.is_none_or(|(b, _)| l >= b) {
                lb = Some((l, i));
            }
        }
        let longest = match lb {
            Some((l, i)) if l >= agg.longest.0 => normal((sb[i].0, sb[i].1)),
            _ => agg.longest.1,
        };
        let width = |d: Dir| -> Option<i128> { Some(known(d)?.2 + ext_of(d, itb()).2) };
        let mut thin = (f64::INFINITY, Dir::X);
        for d in Dir::FIXED.into_iter().chain([longest]) {
            if d.nx == 0 && d.ny == 0 {
                continue;
            }
            let w = width(d)? as f64 / libm::hypot(d.nx as f64, d.ny as f64);
            if w < thin.0 {
                thin = (w, d);
            }
        }
        let thin = thin.1;
        let mut best = ((f64::INFINITY, f64::INFINITY, f64::INFINITY), Dir::X);
        for d in Dir::FIXED.into_iter().chain([thin]) {
            let (alo, ahi, aw) = known(d)?;
            let (blo, bhi, bw) = ext_of(d, itb());
            let inner = (ahi.min(bhi) - alo.max(blo)) as f64;
            let union = (ahi.max(bhi) - alo.min(blo)).max(1) as f64;
            let norm = libm::hypot(d.nx as f64, d.ny as f64);
            let score = (
                inner.max(0.0) / union,
                (inner.min(0.0)) / norm,
                (aw + bw) as f64 / norm,
            );
            if score < best.0 {
                best = (score, d);
            }
        }
        Some(best.1)
    }

    /// Exact minimum distance between the prepared geometry and `b`, with a pair of
    /// closest points: [`distance`](crate::distance)`(g, b)`, the same pair included.
    pub fn distance<B: Geometry + ?Sized>(&self, b: &B) -> Option<Closest> {
        let ba = self.bbox?;
        let bb = b.bbox()?;
        if !self.in_range || !in_range(b) {
            return None;
        }
        let everything = Rect {
            min: Point::new(i64::MIN, i64::MIN),
            max: Point::new(i64::MAX, i64::MAX),
        };
        let mut sb = collect_segments(b, &everything);
        if self.segs.len() + sb.len() < SMALL {
            return crate::distance::distance(self.geom, b);
        }
        if ba.intersects(&bb) {
            // The full search sweeps along `dir` and reports the first touching pair it
            // meets; pairs with meeting boxes involve segments near `b` only, and their
            // order in the sweep does not depend on the other segments.
            let Some(dir) = self.separating(&sb) else {
                return crate::distance::distance(self.geom, b);
            };
            let sbox = sb.iter().fold(bb, |r, s| r.union(&s.2));
            let sa = self.collect(&sbox);
            let mut hit: Option<PointF> = None;
            dir_sweep(&sa, &sb, dir, 0, |x, y| {
                if segments_intersect(x.0, x.1, y.0, y.1) {
                    hit = Some(crate::distance::common_point(x.0, x.1, y.0, y.1));
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
            if let Some(p) = self.component_in(b).or_else(|| self.first_inside(b)) {
                return Some(Closest {
                    sq: SqDist::ZERO,
                    a: p.into(),
                    b: p.into(),
                });
            }
        }
        let bvh = Bvh::build(&mut sb);
        let radius =
            |s: &SqDist| -> u128 { (libm::ceil(s.distance_f64()) as u128).saturating_add(2) };
        let sa0 = self.segs[0];
        let init = segment_segment(sa0.0, sa0.1, sb[0].0, sb[0].1);
        // The minimum, and the first segment (in scanning order) reaching it.
        let mut best = (init.0, 0u32);
        let mut r = radius(&best.0);
        for y in &sb {
            let mut stack = vec![0u32];
            while let Some(ni) = stack.pop() {
                let (nb, first, count) = self.index.nodes[ni as usize];
                if rect_gap2(&y.2, &nb) > r * r {
                    continue;
                }
                if count == 0 {
                    let gl = rect_gap2(&y.2, &self.index.nodes[first as usize].0);
                    let gr = rect_gap2(&y.2, &self.index.nodes[first as usize + 1].0);
                    if gl <= gr {
                        stack.push(first + 1);
                        stack.push(first);
                    } else {
                        stack.push(first);
                        stack.push(first + 1);
                    }
                    continue;
                }
                for &i in &self.index.items[first as usize..(first + count) as usize] {
                    let x = &self.segs[i as usize];
                    if rect_gap2(&x.2, &y.2) > r * r {
                        continue;
                    }
                    let s = segment_segment(x.0, x.1, y.0, y.1).0;
                    match s.cmp(&best.0) {
                        core::cmp::Ordering::Less => {
                            best = (s, i);
                            r = radius(&best.0);
                        }
                        core::cmp::Ordering::Equal if i < best.1 => best.1 = i,
                        _ => {}
                    }
                }
            }
        }
        let min = best.0;
        if init.0 == min {
            return Some(Closest {
                sq: init.0,
                a: init.1,
                b: init.2,
            });
        }
        // The full search keeps the first pair reaching the minimum: the first such
        // segment, against the first segment of `b` reaching it in its traversal order.
        let x = self.segs[best.1 as usize];
        let mut stack: Vec<u32> = vec![0];
        while let Some(ni) = stack.pop() {
            let node = &bvh.nodes[ni as usize];
            if rect_gap2(&x.2, &node.bbox) > r * r {
                continue;
            }
            if node.count > 0 {
                for y in &sb[node.first as usize..(node.first + node.count) as usize] {
                    let (s, p, q) = segment_segment(x.0, x.1, y.0, y.1);
                    if s == min {
                        return Some(Closest { sq: s, a: p, b: q });
                    }
                }
            } else {
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
        // Not reached: the minimum is attained by `x` and some segment of `b`.
        crate::distance::distance(self.geom, b)
    }
}
