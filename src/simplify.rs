//! Tolerance-based simplification (Douglas–Peucker) of paths and polygon sets.
//!
//! [`simplify_path`] is the plain Douglas–Peucker algorithm. [`simplify_polygons`] runs the
//! same recursion on every ring of a polygon set but accepts a shortcut only when replacing
//! the chain it removes cannot change the topology of the set:
//!
//! * the shortcut meets no other current edge, except at its own endpoints without
//!   overlapping (so no new crossing or touching point appears);
//! * the cyclic order of the edges around a shortcut endpoint shared with other rings is
//!   unchanged (a pinched hole stays on the same side);
//! * the region swept between the removed chain and the shortcut contains no other ring
//!   and not the rest of the ring itself (nesting is preserved and a ring never turns inside
//!   out);
//! * the ring keeps at least three vertices and the sign of its area.
//!
//! Each accepted step keeps a valid set valid, so the whole run does. Vertices shared by
//! several rings, or touching another ring's edge, are never removed.

use crate::geom::{Path, Point, Polygon, PolygonSet, Rect, Ring};
use crate::predicates::{
    cmp_angle, cross, dist2, dot, in_segment_interior, on_segment, orient, segments_intersect, sub,
};
use crate::query::{ring_area2, ring_winding};
use core::cmp::Ordering;

/// Tolerances above this behave identically (it exceeds the diameter of the coordinate
/// range); clamping keeps every squared quantity within `u128`.
const MAX_TOL: i64 = 1 << 42;

/// Simplifies an open path with the Douglas–Peucker algorithm.
///
/// The first and last vertices are always kept and every output vertex is an input vertex.
/// Each removed vertex lies within Euclidean distance `tolerance` of the output edge that
/// replaces it (compared exactly), so the Hausdorff distance between input and output is at
/// most `tolerance`. Consecutive duplicate vertices and vertices lying strictly inside the
/// segment joining their neighbours are always removed; with `tolerance <= 0` nothing else
/// is.
///
/// This is not topology preserving: the output may cross itself where the input did not.
/// Paths with coordinates outside `±`[`MAX_COORD`](crate::MAX_COORD) are returned
/// unchanged.
pub fn simplify_path(path: &Path, tolerance: i64) -> Path {
    if path.iter().any(|p| !p.in_range()) {
        return path.clone();
    }
    let tol = tolerance.clamp(0, MAX_TOL);
    let mut pts: Vec<Point> = path.to_vec();
    pts.dedup();
    let n = pts.len();
    if n <= 2 {
        return Path(pts);
    }
    let mut keep = vec![false; n];
    keep[0] = true;
    keep[n - 1] = true;
    let mut stack = vec![(0usize, n - 1)];
    while let Some((i, j)) = stack.pop() {
        if j < i + 2 {
            continue;
        }
        let scan = scan_chain(&pts, i, j, tol);
        if !scan.within {
            keep[scan.far] = true;
            stack.push((scan.far, j));
            stack.push((i, scan.far));
        }
    }
    let mut out: Vec<Point> = Vec::with_capacity(n);
    for (k, &p) in pts.iter().enumerate() {
        if !keep[k] {
            continue;
        }
        while out.len() >= 2 && in_segment_interior(out[out.len() - 2], p, out[out.len() - 1]) {
            out.pop();
        }
        out.push(p);
    }
    Path(out)
}

/// Simplifies a single polygon, preserving its topology.
///
/// Equivalent to [`simplify_polygons`] on a one-element set; see there for the guarantees.
pub fn simplify_polygon(poly: &Polygon, tolerance: i64) -> Polygon {
    simplify_polygons(core::slice::from_ref(poly), tolerance)
        .pop()
        .unwrap_or_default()
}

/// Simplifies every ring of a polygon set with a topology-preserving Douglas–Peucker
/// algorithm, considering all rings together.
///
/// Guarantees, for a valid input set (see [`validate_set`](crate::validate_set)):
///
/// * every output vertex is an input vertex of the same ring, in the same cyclic order, and
///   every removed vertex lies within Euclidean distance `tolerance` (compared exactly) of
///   the output edge replacing it, so each ring moves by at most `tolerance` (Hausdorff);
/// * the output is valid: no edge crosses or touches another edge except where the input
///   already touched, nesting and orientation of every ring are unchanged, no ring is
///   dropped or collapsed (each keeps at least three vertices and a non-zero area of the
///   same sign);
/// * vertices shared by several rings (pinch points) or lying on another ring's edge are
///   kept;
/// * canonical input stays canonical (see [`check_canonical`](crate::check_canonical)):
///   rings keep starting at their minimum vertex, vertices lying strictly between their
///   neighbours are removed unless shared, and holes and polygons are re-sorted when they
///   were sorted on input.
///
/// With `tolerance <= 0` only duplicate and collinear (straight-through) vertices are
/// removed. The result has one polygon per input polygon and one ring per input ring.
///
/// Invalid input never panics: simplification is then best effort (each step still applies
/// the same local checks). Rings with fewer than three distinct vertices or zero area are
/// returned without duplicate vertices but otherwise unchanged, and sets with coordinates
/// outside `±`[`MAX_COORD`](crate::MAX_COORD) are returned unchanged.
///
/// Runs in about `O(n log n)` for typical inputs, using a uniform grid index.
pub fn simplify_polygons(polys: &[Polygon], tolerance: i64) -> PolygonSet {
    let total: usize = polys.iter().map(|p| p.vertex_count()).sum();
    if total == 0
        || total >= (u32::MAX / 4) as usize
        || polys
            .iter()
            .any(|p| p.rings().any(|r| r.iter().any(|q| !q.in_range())))
    {
        return polys.to_vec();
    }
    // Work in a canonical order (rings rotated to their smallest vertex, holes and
    // polygons sorted) so the greedy simplification does not depend on input order or
    // ring start vertices; results are mapped back to the caller's order.
    let rot = |r: &Ring| -> Ring {
        let mut v = r.0.clone();
        if let Some((i, _)) = v.iter().enumerate().min_by_key(|(_, p)| **p) {
            v.rotate_left(i);
        }
        Ring(v)
    };
    let mut work: Vec<(Polygon, usize, Vec<usize>)> = polys
        .iter()
        .enumerate()
        .map(|(pi, p)| {
            let mut holes: Vec<(Ring, usize)> = p
                .holes
                .iter()
                .enumerate()
                .map(|(hi, h)| (rot(h), hi))
                .collect();
            holes.sort_by(|a, b| a.0.0.cmp(&b.0.0).then(a.1.cmp(&b.1)));
            let order = holes.iter().map(|h| h.1).collect();
            (
                Polygon {
                    outer: rot(&p.outer),
                    holes: holes.into_iter().map(|h| h.0).collect(),
                },
                pi,
                order,
            )
        })
        .collect();
    work.sort_by(|a, b| {
        a.0.outer
            .0
            .cmp(&b.0.outer.0)
            .then_with(|| {
                a.0.holes
                    .iter()
                    .map(|h| &h.0)
                    .cmp(b.0.holes.iter().map(|h| &h.0))
            })
            .then(a.1.cmp(&b.1))
    });
    let canon: Vec<Polygon> = work.iter().map(|w| w.0.clone()).collect();
    let mut s = Simplifier::new(&canon, tolerance.clamp(0, MAX_TOL));
    s.run();
    let done = s.output(&canon);
    let mut out: Vec<Polygon> = vec![Polygon::default(); polys.len()];
    for ((_, pi, order), mut p) in work.into_iter().zip(done) {
        let mut holes = vec![Ring::default(); order.len()];
        for (h, hi) in core::mem::take(&mut p.holes).into_iter().zip(order) {
            holes[hi] = h;
        }
        p.holes = holes;
        out[pi] = p;
    }
    // Canonical input stays canonical: keep its sort orders.
    if polys.windows(2).all(|w| w[0].outer.0 < w[1].outer.0) {
        out.sort_by(|a, b| a.outer.0.cmp(&b.outer.0));
    }
    for (p, orig) in out.iter_mut().zip(polys) {
        if orig.holes.windows(2).all(|w| w[0].0 < w[1].0) {
            p.holes.sort_by(|a, b| a.0.cmp(&b.0));
        }
    }
    out
}

/// Result of scanning the interior vertices of a chain against its shortcut.
struct Scan {
    /// Every interior vertex is within tolerance.
    within: bool,
    /// Index of the (approximately) farthest interior vertex.
    far: usize,
    /// Bounding box of the chain (endpoints included).
    bbox: Rect,
}

/// Scans `pts[i+1..j]` against the segment `pts[i]-pts[j]` (`j >= i + 2`).
fn scan_chain(pts: &[Point], i: usize, j: usize, tol: i64) -> Scan {
    let (a, b) = (pts[i], pts[j]);
    let t2 = (tol as u128) * (tol as u128);
    let mut bbox = Rect::new(a, b);
    let mut within = true;
    let mut far = i + 1;
    let mut fd = -1.0f64;
    for (k, &p) in pts.iter().enumerate().take(j).skip(i + 1) {
        bbox.add_point(p);
        let (d, ok) = seg_dist(p, a, b, t2, within);
        if d > fd {
            fd = d;
            far = k;
        }
        within &= ok;
    }
    Scan { within, far, bbox }
}

/// Approximate squared distance from `p` to segment `a-b` and, when `exact` is set, whether
/// that distance is at most `sqrt(t2)` (decided exactly). When `exact` is unset the flag is
/// `false`.
fn seg_dist(p: Point, a: Point, b: Point, t2: u128, exact: bool) -> (f64, bool) {
    if a == b || dot(sub(p, a), sub(b, a)) <= 0 {
        let d = dist2(a, p) as u128;
        return (d as f64, exact && d <= t2);
    }
    if dot(sub(p, b), sub(a, b)) <= 0 {
        let d = dist2(b, p) as u128;
        return (d as f64, exact && d <= t2);
    }
    // Perpendicular distance: |cross| / |ab|; compare cross^2 <= t2 * |ab|^2.
    let c = orient(a, b, p).unsigned_abs();
    let l2 = dist2(a, b) as u128;
    let cf = c as f64;
    let lf = l2 as f64;
    let d = cf * cf / lf;
    if !exact {
        return (d, false);
    }
    let lhs = cf * cf;
    let rhs = t2 as f64 * lf;
    let ok = if lhs < rhs * (1.0 - 1e-9) {
        true
    } else if lhs > rhs * (1.0 + 1e-9) {
        false
    } else {
        mul_wide(c, c) <= mul_wide(t2, l2)
    };
    (d, ok)
}

/// Full 256-bit product of two `u128`, as `(high, low)` (lexicographic order = numeric).
fn mul_wide(a: u128, b: u128) -> (u128, u128) {
    const M: u128 = u64::MAX as u128;
    let (a1, a0) = (a >> 64, a & M);
    let (b1, b0) = (b >> 64, b & M);
    let p00 = a0 * b0;
    let p01 = a0 * b1;
    let p10 = a1 * b0;
    let p11 = a1 * b1;
    let mid = (p00 >> 64) + (p01 & M) + (p10 & M);
    let lo = (p00 & M) | (mid << 64);
    let hi = p11 + (p01 >> 64) + (p10 >> 64) + (mid >> 64);
    (hi, lo)
}

/// `true` when the shortcut `a-b` meets edge `c-d` anywhere other than at a single point that
/// is `a` or `b`.
fn shortcut_conflicts(a: Point, b: Point, c: Point, d: Point) -> bool {
    if !segments_intersect(a, b, c, d) {
        return false;
    }
    if orient(a, b, c) != 0 || orient(a, b, d) != 0 {
        // Not collinear: a single common point, allowed only at a shortcut endpoint.
        return !(on_segment(c, d, a) || on_segment(c, d, b));
    }
    // Collinear: the common part is an interval; allowed only when it is a single point
    // (necessarily `a` or `b`).
    let horiz = (b.x - a.x).abs() >= (b.y - a.y).abs();
    let key = |p: Point| if horiz { p.x } else { p.y };
    let (s0, s1) = (key(a).min(key(b)), key(a).max(key(b)));
    let (t0, t1) = (key(c).min(key(d)), key(c).max(key(d)));
    s1.min(t1) > s0.max(t0)
}

/// Compares directions `v` and `w` by counter-clockwise angle measured from direction `p`.
fn ccw_from(p: Point, v: Point, w: Point) -> Ordering {
    let wv = cmp_angle(v, p) == Ordering::Less;
    let ww = cmp_angle(w, p) == Ordering::Less;
    wv.cmp(&ww).then_with(|| cmp_angle(v, w))
}

/// `true` when direction `d` lies strictly inside the counter-clockwise sweep from `p` to
/// `x`.
fn strictly_between(p: Point, x: Point, d: Point) -> bool {
    ccw_from(p, p, d) == Ordering::Less && ccw_from(p, d, x) == Ordering::Less
}

/// Uniform grid over a bounding box; each cell lists ids: those given when the grid was
/// filled (compactly, in `items[start[c]..start[c + 1]]`), then those added later.
struct Grid {
    x0: i64,
    y0: i64,
    s: i64,
    nx: i64,
    ny: i64,
    start: Vec<u32>,
    items: Vec<u32>,
    extra: Vec<Vec<u32>>,
}

impl Grid {
    /// A grid with about `n` cells over `bb`.
    fn new(bb: Rect, n: usize) -> Grid {
        let w = bb.width() as f64 + 1.0;
        let h = bb.height() as f64 + 1.0;
        let n = n.max(1) as f64;
        let s = (w * h / n).sqrt().max(w.max(h) / n).ceil().max(1.0);
        let s = if s >= (1u64 << 62) as f64 {
            1i64 << 62
        } else {
            s as i64
        };
        let nx = bb.width() / s + 1;
        let ny = bb.height() / s + 1;
        let nc = (nx * ny) as usize;
        Grid {
            x0: bb.min.x,
            y0: bb.min.y,
            s,
            nx,
            ny,
            start: vec![0; nc + 1],
            items: Vec::new(),
            extra: vec![Vec::new(); nc],
        }
    }

    /// Fills the (empty) grid with `(cell, id)` pairs, keeping their order within a cell.
    fn fill(&mut self, pairs: &[(u32, u32)]) {
        for &(c, _) in pairs {
            self.start[c as usize + 1] += 1;
        }
        for c in 0..self.extra.len() {
            self.start[c + 1] += self.start[c];
        }
        let mut pos = self.start.clone();
        self.items = vec![0; pairs.len()];
        for &(c, id) in pairs {
            self.items[pos[c as usize] as usize] = id;
            pos[c as usize] += 1;
        }
    }

    /// The ids listed in cell `c`.
    #[inline]
    fn cell(&self, c: usize) -> impl Iterator<Item = u32> + '_ {
        self.items[self.start[c] as usize..self.start[c + 1] as usize]
            .iter()
            .chain(self.extra[c].iter())
            .copied()
    }

    #[inline]
    fn col(&self, x: i64) -> i64 {
        (x - self.x0).div_euclid(self.s).clamp(0, self.nx - 1)
    }

    #[inline]
    fn row(&self, y: i64) -> i64 {
        (y - self.y0).div_euclid(self.s).clamp(0, self.ny - 1)
    }

    #[inline]
    fn colf(&self, x: f64) -> i64 {
        let c = ((x - self.x0 as f64) / self.s as f64).floor();
        if c < 0.0 {
            0
        } else if c >= self.nx as f64 {
            self.nx - 1
        } else {
            c as i64
        }
    }

    #[inline]
    fn cell_of(&self, p: Point) -> usize {
        (self.row(p.y) * self.nx + self.col(p.x)) as usize
    }

    /// Pushes (a superset of) the cells met by segment `a-b` into `out`.
    fn segment_cells(&self, a: Point, b: Point, out: &mut Vec<usize>) {
        out.clear();
        let (ylo, yhi) = (a.y.min(b.y), a.y.max(b.y));
        let (r0, r1) = (self.row(ylo), self.row(yhi));
        let dy = (b.y - a.y) as f64;
        let dx = (b.x - a.x) as f64;
        for r in r0..=r1 {
            let (c0, c1) = if a.y == b.y || r0 == r1 {
                (self.col(a.x.min(b.x)), self.col(a.x.max(b.x)))
            } else {
                let ry0 = (self.y0 + r * self.s).max(ylo) as f64;
                let ry1 = (self.y0 + (r + 1) * self.s).min(yhi) as f64;
                let x_at = |y: f64| a.x as f64 + dx * ((y - a.y as f64) / dy);
                let (xa, xb) = (x_at(ry0), x_at(ry1));
                (
                    self.colf(xa.min(xb)).saturating_sub(1).max(0),
                    (self.colf(xa.max(xb)) + 1).min(self.nx - 1),
                )
            };
            for c in c0..=c1 {
                out.push((r * self.nx + c) as usize);
            }
        }
    }

    /// Cell index ranges `(c0, c1, r0, r1)` covering `bb`.
    fn rect_range(&self, bb: &Rect) -> (i64, i64, i64, i64) {
        (
            self.col(bb.min.x),
            self.col(bb.max.x),
            self.row(bb.min.y),
            self.row(bb.max.y),
        )
    }
}

/// One edge of the current configuration.
struct Edge {
    a: Point,
    b: Point,
    alive: bool,
}

/// One ring being simplified.
struct RingData {
    /// Deduplicated vertices with the first repeated at the end (`n + 1` entries, or none).
    pts: Vec<Point>,
    n: usize,
    /// Global id of vertex 0 (and of the original edge leaving it).
    base: usize,
    next: Vec<u32>,
    prev: Vec<u32>,
    /// Current vertex count.
    count: usize,
    /// Current twice-signed area.
    area: i128,
    /// Valid enough to simplify (at least three vertices, non-zero area).
    active: bool,
}

struct Simplifier {
    tol: i64,
    rings: Vec<RingData>,
    edges: Vec<Edge>,
    grid: Grid,
    /// Ring ids indexed by the position of their vertex 0 (never removed).
    anchors: Grid,
    pinned: Vec<bool>,
    stamps: Vec<u32>,
    stamp: u32,
    buf: Vec<usize>,
}

impl Simplifier {
    fn new(polys: &[Polygon], tol: i64) -> Simplifier {
        let mut rings = Vec::new();
        let mut base = 0;
        let mut bb: Option<Rect> = None;
        for poly in polys {
            for ring in poly.rings() {
                let mut v = ring.to_vec();
                v.dedup();
                while v.len() > 1 && v.first() == v.last() {
                    v.pop();
                }
                let n = v.len();
                let area = ring_area2(&v);
                for &p in &v {
                    match &mut bb {
                        Some(r) => r.add_point(p),
                        None => bb = Some(Rect::new(p, p)),
                    }
                }
                if n > 0 {
                    v.push(v[0]);
                }
                rings.push(RingData {
                    pts: v,
                    n,
                    base,
                    next: (0..n).map(|k| ((k + 1) % n) as u32).collect(),
                    prev: (0..n).map(|k| ((k + n - 1) % n) as u32).collect(),
                    count: n,
                    area,
                    active: n >= 3 && area != 0,
                });
                base += n;
            }
        }
        let bb = bb.unwrap_or(Rect::new(Point::default(), Point::default()));
        let mut grid = Grid::new(bb, base);
        let mut anchors = Grid::new(bb, base);
        let mut edges = Vec::with_capacity(base + base / 4);
        let mut buf = Vec::new();
        let mut anchor_pairs: Vec<(u32, u32)> = Vec::with_capacity(rings.len());
        let mut pairs: Vec<(u32, u32)> = Vec::with_capacity(base + base / 2);
        for (ri, r) in rings.iter().enumerate() {
            if r.n > 0 {
                anchor_pairs.push((anchors.cell_of(r.pts[0]) as u32, ri as u32));
            }
            for k in 0..r.n {
                let (a, b) = (r.pts[k], r.pts[k + 1]);
                let id = edges.len() as u32;
                edges.push(Edge { a, b, alive: true });
                grid.segment_cells(a, b, &mut buf);
                pairs.extend(buf.iter().map(|&c| (c as u32, id)));
            }
        }
        anchors.fill(&anchor_pairs);
        grid.fill(&pairs);
        drop(pairs);
        // Pinned vertices: shared with another vertex, lying on another edge, or ends of an
        // edge another vertex lies on.
        let mut pinned = vec![false; base];
        let mut all: Vec<(Point, usize)> = Vec::with_capacity(base);
        let mut edge_ring: Vec<usize> = Vec::with_capacity(base);
        for (ri, r) in rings.iter().enumerate() {
            for k in 0..r.n {
                all.push((r.pts[k], r.base + k));
                edge_ring.push(ri);
            }
        }
        let all = crate::par::bucket_sort_by_x(all, |e| e.0.x, |a, b| a.cmp(b));
        let mut i = 0;
        while i < all.len() {
            let mut j = i + 1;
            while j < all.len() && all[j].0 == all[i].0 {
                j += 1;
            }
            if j - i > 1 {
                for e in &all[i..j] {
                    pinned[e.1] = true;
                }
            }
            i = j;
        }
        for &(q, g) in &all {
            for e in grid.cell(grid.cell_of(q)) {
                let e = e as usize;
                let (a, b) = (edges[e].a, edges[e].b);
                // (The box test first: most edges of the cell are far from `q`.)
                if Rect::new(a, b).contains_point(q) && in_segment_interior(a, b, q) {
                    pinned[g] = true;
                    let r = &rings[edge_ring[e]];
                    pinned[e] = true;
                    pinned[r.base + (e - r.base + 1) % r.n] = true;
                }
            }
        }
        let stamps = vec![0; edges.len()];
        Simplifier {
            tol,
            rings,
            edges,
            grid,
            anchors,
            pinned,
            stamps,
            stamp: 0,
            buf,
        }
    }

    fn run(&mut self) {
        let mut stack: Vec<(usize, usize)> = Vec::new();
        for r in 0..self.rings.len() {
            let ring = &self.rings[r];
            if !ring.active {
                continue;
            }
            let n = ring.n;
            let mut anchors = vec![0, n];
            let o = ring.pts[0];
            let far = (1..n).max_by_key(|&k| dist2(o, ring.pts[k])).unwrap_or(0);
            anchors.push(far);
            anchors.extend((1..n).filter(|&k| self.pinned[ring.base + k]));
            anchors.sort_unstable();
            anchors.dedup();
            for w in anchors.windows(2).rev() {
                stack.push((w[0], w[1]));
            }
            while let Some((i, j)) = stack.pop() {
                if let Some(k) = self.step(r, i, j) {
                    stack.push((k, j));
                    stack.push((i, k));
                }
            }
        }
    }

    fn next_stamp(&mut self) -> u32 {
        if self.stamp == u32::MAX {
            self.stamps.iter_mut().for_each(|s| *s = 0);
            self.stamp = 0;
        }
        self.stamp += 1;
        self.stamp
    }

    /// Tries to replace the chain `i..=j` of ring `r` by a single edge. Returns the vertex to
    /// split at when the shortcut is rejected.
    fn step(&mut self, r: usize, i: usize, j: usize) -> Option<usize> {
        if j < i + 2 {
            return None;
        }
        let st = self.next_stamp();
        let tol = self.tol;
        let Simplifier {
            rings,
            edges,
            grid,
            anchors,
            pinned,
            stamps,
            buf,
            ..
        } = self;
        let ring = &rings[r];
        let n = ring.n;
        let (a, b) = (ring.pts[i], ring.pts[j]);
        let scan = scan_chain(&ring.pts, i, j, tol);
        let split = Some(scan.far);
        let removed = j - i - 1;
        if !scan.within || a == b || ring.count < removed + 3 {
            return split;
        }
        let chain = &ring.pts[i..=j];
        let new_area = ring.area - ring_area2(chain);
        if new_area == 0 || (new_area > 0) != (ring.area > 0) {
            return split;
        }
        let jj = j % n;
        // The rest of the ring must stay outside the swept region.
        let q = ring.pts[ring.next[jj] as usize];
        if q == a || q == b || ring_winding(chain, q) != Some(0) {
            return split;
        }
        // No other edge may meet the shortcut (the chain's own edges are being removed).
        let (lo, hi) = (ring.base + i, ring.base + j);
        grid.segment_cells(a, b, buf);
        for &c in buf.iter() {
            for e in grid.cell(c) {
                let e = e as usize;
                if stamps[e] == st {
                    continue;
                }
                stamps[e] = st;
                let ed = &edges[e];
                if !ed.alive || (lo..hi).contains(&e) {
                    continue;
                }
                if shortcut_conflicts(a, b, ed.a, ed.b) {
                    return split;
                }
            }
        }
        // Cyclic order of edges around shared endpoints must be preserved.
        let ends = [
            (i, ring.pts[ring.prev[i] as usize], ring.pts[i + 1], b),
            (jj, ring.pts[ring.next[jj] as usize], ring.pts[j - 1], a),
        ];
        for (k, other, old, new) in ends {
            if !pinned[ring.base + k] {
                continue;
            }
            let v = ring.pts[k];
            let (p, u, s) = (sub(other, v), sub(old, v), sub(new, v));
            for e in grid.cell(grid.cell_of(v)) {
                let e = e as usize;
                let ed = &edges[e];
                if !ed.alive || (lo..hi).contains(&e) || !on_segment(ed.a, ed.b, v) {
                    continue;
                }
                for w in [ed.a, ed.b] {
                    if w == v {
                        continue;
                    }
                    let d = sub(w, v);
                    if cross(d, p) == 0 && dot(d, p) > 0 {
                        continue;
                    }
                    if strictly_between(p, u, d) != strictly_between(p, s, d) {
                        return split;
                    }
                }
            }
        }
        // No other ring may lie in the swept region. Each other ring is entirely inside or
        // outside it (it meets neither the chain nor the shortcut), so one vertex decides.
        let bb = scan.bbox;
        let (c0, c1, r0, r1) = anchors.rect_range(&bb);
        let ncells = ((c1 - c0 + 1) as u128) * ((r1 - r0 + 1) as u128);
        let check = |r2: usize| -> bool {
            let o = &rings[r2];
            if r2 == r || o.n == 0 || !bb.contains_point(o.pts[0]) {
                return true;
            }
            let mut k = 0usize;
            for _ in 0..o.n {
                let p = o.pts[k];
                if p != a && p != b {
                    return !bb.contains_point(p) || ring_winding(chain, p) == Some(0);
                }
                k = o.next.get(k).map_or(0, |&x| x as usize);
            }
            false
        };
        if ncells > rings.len() as u128 {
            if !(0..rings.len()).all(check) {
                return split;
            }
        } else {
            for row in r0..=r1 {
                for col in c0..=c1 {
                    let mut cell = anchors.cell((row * anchors.nx + col) as usize);
                    if !cell.all(|r2| check(r2 as usize)) {
                        return split;
                    }
                }
            }
        }
        // Accept.
        let ring = &mut rings[r];
        ring.next[i] = jj as u32;
        ring.prev[jj] = i as u32;
        ring.count -= removed;
        ring.area = new_area;
        for e in &mut edges[lo..hi] {
            e.alive = false;
        }
        let id = edges.len() as u32;
        edges.push(Edge { a, b, alive: true });
        stamps.push(0);
        grid.segment_cells(a, b, buf);
        for &c in buf.iter() {
            grid.extra[c].push(id);
        }
        None
    }

    /// Output rings in the order of `polys` (re-sorting is left to the caller).
    fn output(&self, polys: &[Polygon]) -> PolygonSet {
        let mut rings = self.rings.iter().map(|r| self.ring_output(r));
        let mut out: PolygonSet = Vec::with_capacity(polys.len());
        for poly in polys {
            let outer = rings.next().unwrap_or_default();
            let holes: Vec<_> = poly
                .holes
                .iter()
                .map(|_| rings.next().unwrap_or_default())
                .collect();
            out.push(Polygon::new(outer, holes));
        }
        out
    }

    fn ring_output(&self, r: &RingData) -> crate::geom::Ring {
        let n = r.n;
        if n == 0 {
            return crate::geom::Ring::new();
        }
        let src = &r.pts[..n];
        if !r.active {
            return src.iter().copied().collect();
        }
        let starts_at_min = src.iter().min() == Some(&src[0]);
        // Current vertices, then drop straight-through vertices that are not pinned.
        let mut idx: Vec<usize> = Vec::with_capacity(r.count);
        let mut k = 0usize;
        for _ in 0..r.count {
            idx.push(k);
            k = r.next[k] as usize;
        }
        let m = idx.len();
        let mut next: Vec<usize> = (0..m).map(|t| (t + 1) % m).collect();
        let mut prev: Vec<usize> = (0..m).map(|t| (t + m - 1) % m).collect();
        let mut alive = vec![true; m];
        let mut count = m;
        let mut work: Vec<usize> = (0..m).rev().collect();
        while let Some(t) = work.pop() {
            if count <= 3 {
                break;
            }
            if !alive[t] || self.pinned[r.base + idx[t]] {
                continue;
            }
            let (p, q) = (prev[t], next[t]);
            if in_segment_interior(src[idx[p]], src[idx[q]], src[idx[t]]) {
                alive[t] = false;
                next[p] = q;
                prev[q] = p;
                count -= 1;
                work.push(q);
                work.push(p);
            }
        }
        let mut pts: Vec<Point> = (0..m).filter(|&t| alive[t]).map(|t| src[idx[t]]).collect();
        if starts_at_min && let Some(mi) = (0..pts.len()).min_by_key(|&t| pts[t]) {
            pts.rotate_left(mi);
        }
        crate::geom::Ring(pts)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::Ring;
    use crate::validate::{check_canonical, validate_set};

    fn p(x: i64, y: i64) -> Point {
        Point::new(x, y)
    }

    fn ring(v: &[(i64, i64)]) -> Ring {
        v.iter().map(|&(x, y)| p(x, y)).collect()
    }

    /// The compact grid lists the same ids in the same order as one vector per cell (the
    /// previous layout), including ids added after filling.
    #[test]
    fn grid_matches_vector_per_cell() {
        let mut s: u64 = 17;
        let mut rnd = |m: i64| {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((s >> 33) as i64).rem_euclid(m)
        };
        for it in 0..200 {
            let range = [10i64, 1000, 1 << 40][it % 3];
            let segs: Vec<(Point, Point)> = (0..1 + rnd(300))
                .map(|_| (p(rnd(range), rnd(range)), p(rnd(range), rnd(range))))
                .collect();
            let bb = Rect::new(p(0, 0), p(range, range));
            let mut g = Grid::new(bb, segs.len() / 2);
            let mut reference: Vec<Vec<u32>> = vec![Vec::new(); g.extra.len()];
            let (mut buf, mut pairs) = (Vec::new(), Vec::new());
            let first = segs.len() * 2 / 3;
            for (id, &(a, b)) in segs[..first].iter().enumerate() {
                g.segment_cells(a, b, &mut buf);
                for &c in &buf {
                    pairs.push((c as u32, id as u32));
                    reference[c].push(id as u32);
                }
            }
            g.fill(&pairs);
            for (id, &(a, b)) in segs.iter().enumerate().skip(first) {
                g.segment_cells(a, b, &mut buf);
                for &c in &buf {
                    g.extra[c].push(id as u32);
                    reference[c].push(id as u32);
                }
            }
            for (c, want) in reference.iter().enumerate() {
                assert_eq!(&g.cell(c).collect::<Vec<_>>(), want);
            }
        }
    }

    /// A CCW square of side `s` at `(x0, y0)` whose edges carry inward zig-zag noise of
    /// amplitude `amp` every `step` units (the corners stay exact).
    fn noisy_square(x0: i64, y0: i64, s: i64, step: i64, amp: i64) -> Ring {
        let mut v = Vec::new();
        for side in 0..4 {
            let mut t = 0;
            let mut k = 0;
            while t < s {
                let d = if k % 2 == 1 { amp } else { 0 };
                v.push(match side {
                    0 => p(x0 + t, y0 + d),
                    1 => p(x0 + s - d, y0 + t),
                    2 => p(x0 + s - t, y0 + s - d),
                    _ => p(x0 + d, y0 + s - t),
                });
                t += step;
                k += 1;
            }
        }
        Ring(v)
    }

    fn canon(rings: &[Ring]) -> PolygonSet {
        crate::union_all(rings, crate::FillRule::NonZero).unwrap()
    }

    #[test]
    fn wide_mul() {
        assert_eq!(mul_wide(u128::MAX, u128::MAX), (u128::MAX - 1, 1));
        assert_eq!(mul_wide(1 << 100, 1 << 100), (1 << 72, 0));
        assert_eq!(mul_wide(3, 5), (0, 15));
    }

    #[test]
    fn path_basic() {
        let path = Path::from([(0, 0), (5, 1), (10, 0), (10, 0), (20, 0), (20, 10)]);
        assert_eq!(
            simplify_path(&path, 1),
            Path::from([(0, 0), (20, 0), (20, 10)])
        );
        assert_eq!(
            simplify_path(&path, 0),
            Path::from([(0, 0), (5, 1), (10, 0), (20, 0), (20, 10)])
        );
        // Collinear backtracking is not straight-through and stays.
        let spike = Path::from([(0, 0), (10, 0), (5, 0)]);
        assert_eq!(simplify_path(&spike, 0), spike);
        assert_eq!(simplify_path(&Path::new(), 5), Path::new());
        let m = crate::MAX_COORD;
        let big = Path::from([(-m, -m), (0, 1), (m, m)]);
        assert_eq!(simplify_path(&big, i64::MAX).len(), 2);
        assert_eq!(simplify_path(&big, 0).len(), 3);
    }

    #[test]
    fn exact_tolerance_boundary() {
        // (5, 3) is exactly 3 from the segment (0,0)-(10,0).
        let path = Path::from([(0, 0), (5, 3), (10, 0)]);
        assert_eq!(simplify_path(&path, 3).len(), 2);
        assert_eq!(simplify_path(&path, 2).len(), 3);
        // Distance 5 from a diagonal at extreme coordinates: (3,4) offsets.
        let m = crate::MAX_COORD;
        let path = Path::from([(-m + 4, -m + 4), (-4, 1), (m - 4, m - 4)]);
        // Exact distance from (-4, 1) to y = x is 5/sqrt(2) ~ 3.54.
        assert_eq!(simplify_path(&path, 4).len(), 2);
        assert_eq!(simplify_path(&path, 3).len(), 3);
    }

    #[test]
    fn noisy_square_to_square() {
        let set = canon(&[noisy_square(0, 0, 1000, 50, 3)]);
        assert!(set[0].outer.len() > 40);
        let out = simplify_polygons(&set, 5);
        assert_eq!(check_canonical(&out, true), Ok(()));
        assert_eq!(
            out[0].outer,
            ring(&[(0, 0), (1000, 0), (1000, 1000), (0, 1000)])
        );
        // Tolerance too small: nothing beyond collinear removal.
        let out = simplify_polygons(&set, 0);
        assert_eq!(out, set);
    }

    #[test]
    fn collinear_removed_at_zero_tolerance() {
        let poly = Polygon::new(
            ring(&[(0, 0), (5, 0), (10, 0), (10, 5), (10, 10), (0, 10)]),
            vec![],
        );
        let out = simplify_polygon(&poly, 0);
        assert_eq!(out.outer, ring(&[(0, 0), (10, 0), (10, 10), (0, 10)]));
        // Negative tolerance behaves as zero.
        assert_eq!(simplify_polygon(&poly, -7), out);
    }

    #[test]
    fn holes_preserved() {
        let set = canon(&[noisy_square(0, 0, 1000, 50, 3)]);
        let hole = canon(&[noisy_square(400, 400, 200, 20, 2)]);
        let poly =
            crate::boolean(crate::Op::Difference, &set, &hole, crate::FillRule::NonZero).unwrap();
        assert_eq!(poly[0].holes.len(), 1);
        for tol in [0, 1, 3, 10, 100, 1000, 100_000] {
            let out = simplify_polygons(&poly, tol);
            assert_eq!(check_canonical(&out, true), Ok(()), "tol {tol}");
            assert_eq!(out.len(), 1);
            assert_eq!(out[0].holes.len(), 1);
            assert!(out[0].holes[0].len() >= 3);
        }
    }

    #[test]
    fn sliver_not_collapsed() {
        let poly = Polygon::new(ring(&[(0, 0), (1000, 0), (1000, 1), (0, 1)]), vec![]);
        for tol in [1, 10, 1 << 50, i64::MAX] {
            let out = simplify_polygon(&poly, tol);
            assert_eq!(out.outer.len(), 3);
            assert!(out.outer.signed_area2() > 0);
            assert_eq!(validate_set(&[out]), Ok(()));
        }
    }

    #[test]
    fn close_polygons_not_merged() {
        // Two noisy squares separated by a gap of 4 (noise amplitude 3, tolerance 20).
        let set = canon(&[
            noisy_square(0, 0, 1000, 50, 3),
            noisy_square(1010, 0, 1000, 50, 3),
        ]);
        assert_eq!(set.len(), 2);
        for tol in [5, 20, 200, 5000] {
            let out = simplify_polygons(&set, tol);
            assert_eq!(check_canonical(&out, true), Ok(()), "tol {tol}");
            assert_eq!(out.len(), 2);
        }
    }

    #[test]
    fn small_island_not_jumped() {
        // A U-shaped outer ring wrapped around a small separate square: shortcutting the U
        // would swallow the square.
        let u = ring(&[
            (0, 0),
            (100, 0),
            (100, 30),
            (60, 30),
            (60, 10),
            (40, 10),
            (40, 30),
            (0, 30),
        ]);
        let island = ring(&[(45, 20), (55, 20), (55, 25), (45, 25)]);
        let set = canon(&[u, island]);
        assert_eq!(set.len(), 2);
        let out = simplify_polygons(&set, 1000);
        assert_eq!(check_canonical(&out, true), Ok(()));
        assert_eq!(out.len(), 2);
        // The island is still inside the notch of the U.
        let isl = out.iter().find(|q| q.outer.contains(&p(45, 20))).unwrap();
        assert!(isl.outer.len() >= 3);
        let u = out.iter().find(|q| q.outer.contains(&p(0, 0))).unwrap();
        assert_eq!(
            crate::locate_in_ring(&u.outer, p(50, 22)),
            crate::Location::Outside
        );
    }

    #[test]
    fn pinch_vertex_kept() {
        // Hole touching the outer ring at (50, 0).
        let outer = ring(&[(0, 0), (50, 0), (100, 0), (100, 100), (0, 100)]);
        let hole = ring(&[(50, 0), (40, 30), (50, 60), (60, 30)]);
        let poly = Polygon::new(outer, vec![hole]);
        assert_eq!(validate_set(core::slice::from_ref(&poly)), Ok(()));
        for tol in [0, 5, 50, 1000] {
            let out = simplify_polygon(&poly, tol);
            assert_eq!(
                validate_set(core::slice::from_ref(&out)),
                Ok(()),
                "tol {tol}"
            );
            assert!(out.outer.contains(&p(50, 0)));
            assert!(out.holes[0].contains(&p(50, 0)));
            assert!(out.holes[0].len() >= 3);
        }
    }

    #[test]
    fn invalid_input_no_panic() {
        let m = crate::MAX_COORD;
        let bad = vec![
            Polygon::new(ring(&[(0, 0), (10, 10), (10, 0), (0, 10)]), vec![]),
            Polygon::new(
                ring(&[(0, 0), (0, 0), (0, 0)]),
                vec![ring(&[]), ring(&[(1, 1)])],
            ),
            Polygon::new(ring(&[(0, 0), (5, 0), (10, 0)]), vec![]),
            Polygon::new(ring(&[(-m, -m), (m, -m), (m, m), (0, 3), (-m, m)]), vec![]),
            Polygon::new(
                ring(&[(0, 0), (10, 0), (10, 10), (5, 0), (0, 10)]),
                vec![ring(&[(1, 1), (20, 1), (2, 2)])],
            ),
        ];
        for tol in [0, 1, 3, 100, i64::MAX, i64::MIN] {
            let out = simplify_polygons(&bad, tol);
            assert_eq!(out.len(), bad.len());
        }
        let out_of_range = vec![Polygon::new(ring(&[(0, 0), (m + 1, 0), (0, 5)]), vec![])];
        assert_eq!(simplify_polygons(&out_of_range, 10), out_of_range);
        assert!(simplify_polygons(&[], 10).is_empty());
    }
}
