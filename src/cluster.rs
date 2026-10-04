//! Clustered booleans: groups of input rings that cannot interact are computed independently.
//!
//! Snap rounding is local to the bounding boxes of the segments: every hot pixel centre (an
//! endpoint, or a crossing rounded to the nearest integer point, hence within the integer
//! boxes of both segments) lies in the box of the segments producing it, and a segment
//! meeting a pixel, or carrying its centre, has the centre in its own box. Rings are grouped
//! into clusters: two rings are in the same cluster when the boxes of some of their edges
//! come within [`MARGIN`] of each other (transitively). The union `T(C)` of the edge boxes
//! of a cluster `C` is connected and meets no edge of any other cluster. Then:
//!
//! * the noded arrangement of all rings is the union of the arrangements of the clusters
//!   (no hot pixel, crossing or T-junction of one cluster lies in another's boxes), and the
//!   noded edges of `C` stay in `T(C)`;
//! * every ring of another cluster winds a constant number of times around all of `T(C)`:
//!   its contribution to the winding numbers along `C`'s edges (the *base winding* of `C`)
//!   is found by locating one point of `C`;
//! * likewise every output ring of another cluster contains either all of `T(C)` or none of
//!   it, so it never lies between two nested output rings of `C`: nesting within a cluster is
//!   local, and the parent of a cluster's top-level rings is the innermost output ring of the
//!   other clusters containing a point of the cluster.
//!
//! So rings are computed cluster by cluster (in parallel with the `rayon` feature), a lone
//! simple ring directly, and the canonical tree built from all of them is the one the whole
//! input computed at once gives (tests compare both). Cluster detection is cheap (ring
//! boxes, then boxes of chunks of consecutive edges) and gives up early when ring boxes
//! overlap too much for it to pay off.

use crate::arrangement::InEdge;
use crate::assemble::{RawRing, assemble_rings, canonical_tree, remove_collinear, rotate_to_min};
use crate::boolean::compute;
use crate::geom::{Point, PolyTree, Rect};
use crate::predicates::{dot, orient, segments_intersect, sub};
use crate::query::ring_area2;

/// Edges whose boxes come within this distance are in the same cluster. Boxes merely not
/// meeting would do (see above); one unit of margin is kept for good measure.
const MARGIN: i64 = 1;
/// Smaller inputs are always computed in one piece.
const MIN_EDGES: usize = 1024;
/// Consecutive edges per chunk box.
const CHUNK: usize = 16;
/// Edges per block in the simplicity check of lone rings.
const BLOCK: usize = 8;

/// Union-find over ring indices (smallest index as root).
struct UnionFind(Vec<u32>);

impl UnionFind {
    fn find(&mut self, mut x: u32) -> u32 {
        while self.0[x as usize] != x {
            let p = self.0[self.0[x as usize] as usize];
            self.0[x as usize] = p;
            x = p;
        }
        x
    }

    fn union(&mut self, a: u32, b: u32) {
        let (a, b) = (self.find(a), self.find(b));
        if a != b {
            let (lo, hi) = if a < b { (a, b) } else { (b, a) };
            self.0[hi as usize] = lo;
        }
    }
}

/// The input rings: edge ranges, bounding boxes, and boxes of chunks of [`CHUNK`]
/// consecutive edges.
struct Rings<'a> {
    edges: &'a [InEdge],
    start: Vec<u32>,
    bbox: Vec<Rect>,
    chunk_start: Vec<u32>,
    chunks: Vec<Rect>,
}

impl<'a> Rings<'a> {
    fn new(edges: &'a [InEdge], ring_starts: &[u32]) -> Self {
        let mut start = Vec::with_capacity(ring_starts.len() + 1);
        let mut bbox = Vec::with_capacity(ring_starts.len());
        let mut chunk_start = Vec::with_capacity(ring_starts.len() + 1);
        let mut chunks = Vec::with_capacity(edges.len() / CHUNK + ring_starts.len());
        chunk_start.push(0);
        for (i, &s) in ring_starts.iter().enumerate() {
            let e = ring_starts.get(i + 1).map_or(edges.len(), |&e| e as usize);
            let s = s as usize;
            if e == s {
                continue;
            }
            start.push(s as u32);
            let mut rb = Rect::new(edges[s].a, edges[s].b);
            for c in edges[s..e].chunks(CHUNK) {
                let mut b = Rect::new(c[0].a, c[0].b);
                for x in &c[1..] {
                    b.add_point(x.b);
                }
                rb = rb.union(&b);
                chunks.push(b);
            }
            bbox.push(rb);
            chunk_start.push(chunks.len() as u32);
        }
        start.push(edges.len() as u32);
        Rings {
            edges,
            start,
            bbox,
            chunk_start,
            chunks,
        }
    }

    fn len(&self) -> usize {
        self.bbox.len()
    }

    /// Edges of ring `r`.
    fn edges(&self, r: usize) -> &'a [InEdge] {
        &self.edges[self.start[r] as usize..self.start[r + 1] as usize]
    }

    /// Boxes of the edges of ring `r` meeting `g`, grown by `grow`. Returns the work done
    /// (boxes tested).
    fn edge_boxes(&self, r: usize, g: &Rect, grow: i64, out: &mut Vec<Rect>) -> usize {
        out.clear();
        let (s, e) = (
            self.chunk_start[r] as usize,
            self.chunk_start[r + 1] as usize,
        );
        let edges = self.edges(r);
        let mut work = e - s;
        for c in s..e {
            if !self.chunks[c].intersects(g) {
                continue;
            }
            let k = (c - s) * CHUNK;
            let ch = &edges[k..(k + CHUNK).min(edges.len())];
            work += ch.len();
            out.extend(
                ch.iter()
                    .map(|e| Rect::new(e.a, e.b))
                    .filter(|b| b.intersects(g))
                    .map(|b| b.expand(grow)),
            );
        }
        work
    }

    /// Whether some edges of rings `a` and `b` have boxes within [`MARGIN`]. Adds the work
    /// done (boxes tested) to `work`.
    fn near(
        &self,
        a: usize,
        b: usize,
        work: &mut usize,
        la: &mut Vec<Rect>,
        lb: &mut Vec<Rect>,
    ) -> bool {
        let ga = self.bbox[a].expand(MARGIN);
        let gb = self.bbox[b].expand(MARGIN);
        *work += self.edge_boxes(a, &gb, MARGIN, la);
        if la.is_empty() {
            return false;
        }
        *work += self.edge_boxes(b, &ga, 0, lb);
        if lb.is_empty() {
            return false;
        }
        if la.len() * lb.len() <= 64 {
            *work += la.len() * lb.len();
            return la.iter().any(|x| lb.iter().any(|y| x.intersects(y)));
        }
        // Sweep along x over both lists.
        la.sort_unstable_by_key(|b| b.min.x);
        lb.sort_unstable_by_key(|b| b.min.x);
        let (mut i, mut j) = (0, 0);
        while i < la.len() && j < lb.len() {
            // Test the box starting first against the other list's boxes starting before its end.
            let (x, rest) = if la[i].min.x <= lb[j].min.x {
                i += 1;
                (&la[i - 1], &lb[j..])
            } else {
                j += 1;
                (&lb[j - 1], &la[i..])
            };
            for y in rest {
                *work += 1;
                if y.min.x > x.max.x {
                    break;
                }
                if x.intersects(y) {
                    return true;
                }
            }
        }
        false
    }
}

/// Clusters of rings: `of_ring[r]` and the rings of every cluster (CSR, increasing).
struct Clusters {
    of_ring: Vec<u32>,
    start: Vec<u32>,
    rings: Vec<u32>,
}

/// Groups the rings into clusters, or `None` when that does not look worthwhile (ring boxes
/// overlapping too much, too much work, a single cluster).
fn find_clusters(rings: &Rings<'_>, force: bool) -> Option<Clusters> {
    let n = rings.len();
    let total = rings.bbox.iter().copied().reduce(|a, b| a.union(&b))?;
    let area = |b: &Rect| (b.width() + 2 * MARGIN) as f64 * (b.height() + 2 * MARGIN) as f64;
    // Heavily overlapping boxes: one big cluster, or too many candidate pairs.
    if !force && rings.bbox.iter().map(area).sum::<f64>() > 3.0 * area(&total) {
        return None;
    }
    // Sort and sweep along the axis where the boxes are thinner.
    let (sw, sh) = rings.bbox.iter().fold((0f64, 0f64), |(w, h), b| {
        (w + b.width() as f64, h + b.height() as f64)
    });
    let axis = |p: Point| if sw <= sh { p.x } else { p.y };
    let mut order: Vec<u32> = (0..n as u32).collect();
    order.sort_unstable_by_key(|&r| (axis(rings.bbox[r as usize].min), r));
    let budget = 8 * rings.edges.len() + 4096;
    let mut work = 0usize;
    let mut uf = UnionFind((0..n as u32).collect());
    let mut active: Vec<u32> = Vec::new();
    let (mut ca, mut cb) = (Vec::new(), Vec::new());
    for &r in &order {
        let b = rings.bbox[r as usize].expand(MARGIN);
        let lo = axis(b.min);
        active.retain(|&o| axis(rings.bbox[o as usize].max) >= lo);
        work += active.len();
        for &o in &active {
            if !b.intersects(&rings.bbox[o as usize]) || uf.find(o) == uf.find(r) {
                continue;
            }
            if rings.near(r as usize, o as usize, &mut work, &mut ca, &mut cb) {
                uf.union(r, o);
            }
        }
        if work > budget && !force {
            return None;
        }
        active.push(r);
    }
    // Number clusters by their first ring.
    let mut id = vec![u32::MAX; n];
    let mut of_ring = vec![0u32; n];
    let mut count = 0u32;
    for r in 0..n as u32 {
        let root = uf.find(r) as usize;
        if id[root] == u32::MAX {
            id[root] = count;
            count += 1;
        }
        of_ring[r as usize] = id[root];
    }
    if count < 2 && !force {
        return None;
    }
    let mut size = vec![0usize; count as usize];
    let mut start = vec![0u32; count as usize + 1];
    for (r, &c) in of_ring.iter().enumerate() {
        size[c as usize] += rings.edges(r).len();
        start[c as usize + 1] += 1;
    }
    // One cluster with almost everything: not worth the bookkeeping.
    if !force && size.iter().max().copied().unwrap_or(0) * 10 > rings.edges.len() * 9 {
        return None;
    }
    for c in 0..count as usize {
        start[c + 1] += start[c];
    }
    let mut pos = start.clone();
    let mut list = vec![0u32; n];
    for (r, &c) in of_ring.iter().enumerate() {
        list[pos[c as usize] as usize] = r as u32;
        pos[c as usize] += 1;
    }
    Some(Clusters {
        of_ring,
        start,
        rings: list,
    })
}

/// A static k-d tree over points tagged with a cluster, for box queries. Every node splits
/// its points at the median along the axis where they spread most; `y[k]` tells the axis
/// of the node whose median is at position `k`.
struct KdTree {
    pts: Vec<(Point, u32)>,
    y: Vec<bool>,
}

const KD_LEAF: usize = 8;

impl KdTree {
    fn new(mut pts: Vec<(Point, u32)>) -> Self {
        let mut y = vec![false; pts.len()];
        Self::build(&mut pts, &mut y);
        KdTree { pts, y }
    }

    #[inline]
    fn key(p: Point, y: bool) -> i64 {
        if y { p.y } else { p.x }
    }

    fn build(s: &mut [(Point, u32)], ys: &mut [bool]) {
        if s.len() <= KD_LEAF {
            return;
        }
        let b = Rect::of_points(s.iter().map(|q| &q.0)).unwrap();
        let y = b.height() > b.width();
        let mid = s.len() / 2;
        s.select_nth_unstable_by_key(mid, |q| Self::key(q.0, y));
        ys[mid] = y;
        let (l, r) = s.split_at_mut(mid);
        let (yl, yr) = ys.split_at_mut(mid);
        Self::build(l, yl);
        Self::build(&mut r[1..], &mut yr[1..]);
    }

    /// Calls `f` for every point in the closed box `r`.
    fn query(&self, r: &Rect, f: &mut impl FnMut(Point, u32)) {
        Self::query_in(&self.pts, &self.y, r, f)
    }

    fn query_in(s: &[(Point, u32)], ys: &[bool], r: &Rect, f: &mut impl FnMut(Point, u32)) {
        if s.len() <= KD_LEAF {
            for &(p, c) in s {
                if r.contains_point(p) {
                    f(p, c);
                }
            }
            return;
        }
        let mid = s.len() / 2;
        let (p, c) = s[mid];
        let y = ys[mid];
        let k = Self::key(p, y);
        if Self::key(r.min, y) <= k {
            Self::query_in(&s[..mid], &ys[..mid], r, f);
        }
        if r.contains_point(p) {
            f(p, c);
        }
        if k <= Self::key(r.max, y) {
            Self::query_in(&s[mid + 1..], &ys[mid + 1..], r, f);
        }
    }
}

/// Winding numbers of the ring with edges `edge(0..m)` around the points of `hits` (none
/// on the ring), reported as `f(cluster, winding)` for non-zero ones. Reorders `hits`.
fn windings(
    m: usize,
    edge: impl Fn(usize) -> (Point, Point),
    hits: &mut [(Point, u32)],
    mut f: impl FnMut(u32, i32),
) {
    let side = |a: Point, b: Point, q: Point| -> i32 {
        if a.y <= q.y {
            (b.y > q.y && orient(a, b, q) > 0) as i32
        } else {
            -((b.y <= q.y && orient(a, b, q) < 0) as i32)
        }
    };
    if hits.len() <= 4 {
        for &(q, c) in hits.iter() {
            let w: i32 = (0..m).map(&edge).map(|(a, b)| side(a, b, q)).sum();
            if w != 0 {
                f(c, w);
            }
        }
        return;
    }
    // Many points: only edges spanning a point's height can count for it.
    hits.sort_unstable_by_key(|h| (h.0.y, h.0.x));
    let mut w = vec![0i32; hits.len()];
    for i in 0..m {
        let (a, b) = edge(i);
        let (lo, hi) = (a.y.min(b.y), a.y.max(b.y));
        if lo == hi {
            continue;
        }
        let from = hits.partition_point(|h| h.0.y < lo);
        let to = from + hits[from..].partition_point(|h| h.0.y < hi);
        for k in from..to {
            w[k] += side(a, b, hits[k].0);
        }
    }
    for (k, &(_, c)) in hits.iter().enumerate() {
        if w[k] != 0 {
            f(c, w[k]);
        }
    }
}

/// `true` when the closed polyline through `p` (consecutive points distinct) is a simple
/// ring: non-adjacent edges never meet and adjacent ones only at their common vertex.
///
/// The ring is split into lexicographically monotone chains, whose own edges cannot meet.
/// Pairs of chains with overlapping ranges are compared by a merge along their common range,
/// first by blocks of [`BLOCK`] edges, then edge by edge where block boxes overlap. Also
/// `false` when that would take too long.
fn ring_is_simple(p: &[Point]) -> bool {
    let m = p.len();
    if m < 3 {
        return false;
    }
    let wrap = |i: usize| if i >= m { i - m } else { i };
    let up = |i: usize| p[i] < p[wrap(i + 1)];
    let Some(s) = (0..m).find(|&i| up(i) != up(wrap(i + m - 1))) else {
        return false;
    };
    // Chains in increasing order: vertices (flat), the ring index of the edge from each
    // vertex to the next, and boxes of blocks of edges.
    struct Chain {
        off: usize,
        len: usize,
        boff: usize,
        bbox: Rect,
    }
    let mut chains: Vec<Chain> = Vec::new();
    let mut vs: Vec<Point> = Vec::with_capacity(m + 64);
    let mut ids: Vec<u32> = Vec::with_capacity(m + 64);
    let mut blocks: Vec<Rect> = Vec::with_capacity(m / BLOCK + 64);
    let (mut i, mut done) = (s, 0);
    while done < m {
        let dir = up(i);
        let off = vs.len();
        loop {
            vs.push(p[i]);
            ids.push(i as u32);
            i = wrap(i + 1);
            done += 1;
            if done == m || up(i) != dir {
                break;
            }
        }
        vs.push(p[i]);
        ids.push(u32::MAX);
        let len = vs.len() - off - 1;
        if !dir {
            vs[off..].reverse();
            // The edge between vertices t and t + 1 is now the one listed at t + 1.
            ids[off..].reverse();
            ids.copy_within(off + 1.., off);
            *ids.last_mut().unwrap() = u32::MAX;
        }
        let boff = blocks.len();
        let mut bbox = Rect::new(vs[off], vs[off]);
        for b in (0..len).step_by(BLOCK) {
            let r = Rect::of_points(&vs[off + b..=off + (b + BLOCK).min(len)]).unwrap();
            bbox = bbox.union(&r);
            blocks.push(r);
        }
        chains.push(Chain {
            off,
            len,
            boff,
            bbox,
        });
    }
    // Adjacent edges may share their common vertex, unless they fold back onto each other.
    let adjacent_ok = |e: u32, f: u32| {
        let (e, f) = (e as usize, f as usize);
        let (e, f) = if wrap(e + 1) == f {
            (e, f)
        } else if wrap(f + 1) == e {
            (f, e)
        } else {
            return false;
        };
        let (u, v, w) = (p[e], p[f], p[wrap(f + 1)]);
        orient(u, v, w) != 0 || dot(sub(u, v), sub(w, v)) < 0
    };
    let mut budget = 32 * m + 1024;
    // Merges edges `ia` of chain `a` with edges `ib` of chain `b` (ranges of edge positions),
    // calling `f` on every pair whose ranges overlap (and a few more) until it returns false.
    fn merge(
        vs: &[Point],
        a: &Chain,
        ia: (usize, usize),
        b: &Chain,
        ib: (usize, usize),
        step: usize,
        mut f: impl FnMut(usize, usize) -> bool,
    ) -> bool {
        let at = |c: &Chain, t: usize| vs[c.off + (t * step).min(c.len)];
        let (mut i, mut j) = (ia.0, ib.0);
        while i < ia.1 && j < ib.1 {
            if !f(i, j) {
                return false;
            }
            let (ea, eb) = (at(a, i + 1), at(b, j + 1));
            if ea <= eb {
                i += 1;
            }
            if eb <= ea {
                j += 1;
            }
        }
        true
    }
    let test_edges = |ca: &Chain, i: usize, cb: &Chain, j: usize| -> bool {
        let (a0, a1) = (vs[ca.off + i], vs[ca.off + i + 1]);
        let (b0, b1) = (vs[cb.off + j], vs[cb.off + j + 1]);
        !(segments_intersect(a0, a1, b0, b1) && !adjacent_ok(ids[ca.off + i], ids[cb.off + j]))
    };
    // Pairs of chains with overlapping ranges, by a sweep over their lowest vertices.
    chains.sort_unstable_by_key(|c| vs[c.off]);
    let highest = |c: &Chain| vs[c.off + c.len];
    let mut active: Vec<usize> = Vec::new();
    for y in 0..chains.len() {
        let cb = &chains[y];
        let lo = vs[cb.off];
        active.retain(|&x| highest(&chains[x]) >= lo);
        budget = match budget.checked_sub(active.len()) {
            Some(b) => b,
            None => return false,
        };
        for &x in &active {
            let ca = &chains[x];
            if !ca.bbox.intersects(&cb.bbox) {
                continue;
            }
            // First block of each chain reaching `lo`.
            let first = |c: &Chain| {
                let nb = c.len.div_ceil(BLOCK);
                let ends = |t: usize| vs[c.off + ((t + 1) * BLOCK).min(c.len)];
                let (mut a, mut b) = (0, nb - 1);
                while a < b {
                    let mid = (a + b) / 2;
                    if ends(mid) < lo { a = mid + 1 } else { b = mid }
                }
                a
            };
            let (na, nb) = (ca.len.div_ceil(BLOCK), cb.len.div_ceil(BLOCK));
            let ok = merge(
                &vs,
                ca,
                (first(ca), na),
                cb,
                (first(cb), nb),
                BLOCK,
                |bi, bj| {
                    if budget == 0 {
                        return false;
                    }
                    budget -= 1;
                    if !blocks[ca.boff + bi].intersects(&blocks[cb.boff + bj]) {
                        return true;
                    }
                    let ra = (bi * BLOCK, ((bi + 1) * BLOCK).min(ca.len));
                    let rb = (bj * BLOCK, ((bj + 1) * BLOCK).min(cb.len));
                    merge(&vs, ca, ra, cb, rb, 1, |i, j| {
                        if budget == 0 {
                            return false;
                        }
                        budget -= 1;
                        test_edges(ca, i, cb, j)
                    })
                },
            );
            if !ok {
                return false;
            }
        }
        active.push(y);
    }
    true
}

/// The rings of one cluster: canonical rings, hole flags, parents within the cluster.
#[derive(Default)]
struct Part {
    rings: Vec<RawRing>,
    is_hole: Vec<bool>,
    parent: Vec<Option<u32>>,
}

/// The result of a lone simple ring with base winding `base`, as the full pipeline would
/// give it (`None` when the ring is not simple).
fn lone_ring(
    edges: &[InEdge],
    base: [i32; 2],
    inside: &impl Fn([i32; 2]) -> bool,
    keep_collinear: bool,
) -> Option<Part> {
    let m = edges.len();
    let pts: Vec<Point> = edges.iter().map(|e| e.a).collect();
    if !ring_is_simple(&pts) {
        return None;
    }
    // Collinear vertices; the orientation of a simple ring is the turn at its smallest
    // vertex (a convex hull vertex, never collinear with its neighbours).
    let mut collinear = false;
    let (mut prev, mut low) = (pts[m - 1], 0);
    for (i, &v) in pts.iter().enumerate() {
        let next = pts[if i + 1 == m { 0 } else { i + 1 }];
        collinear |= orient(prev, v, next) == 0;
        if v < pts[low] {
            low = i;
        }
        prev = v;
    }
    let ccw = orient(pts[(low + m - 1) % m], pts[low], pts[(low + 1) % m]) > 0;
    let mut w_in = base;
    w_in[edges[0].operand as usize] += if ccw { 1 } else { -1 };
    let (i_in, i_out) = (inside(w_in), inside(base));
    if i_in == i_out {
        return Some(Part::default());
    }
    // Interior on the left: counter-clockwise when the material is inside.
    let mut r = if ccw == i_in {
        RawRing {
            pts,
            tags: edges.iter().map(|e| e.tag).collect(),
        }
    } else {
        RawRing {
            pts: (0..m).map(|k| pts[(m - k) % m]).collect(),
            tags: (0..m).map(|k| edges[m - 1 - k].tag).collect(),
        }
    };
    // Without collinear vertices there is nothing to remove.
    if !keep_collinear && collinear {
        remove_collinear(&mut r, &[]);
    }
    rotate_to_min(&mut r);
    Some(Part {
        rings: vec![r],
        is_hole: vec![!i_in],
        parent: vec![None],
    })
}

/// Computes the boolean of `edges` (rings starting at `ring_starts`) cluster by cluster, or
/// returns `None` when the input does not split into several clusters (or is too small to
/// bother). The result is identical to the one-piece computation.
pub(crate) fn execute(
    edges: &[InEdge],
    ring_starts: &[u32],
    inside: impl Fn([i32; 2]) -> bool + Sync + Send,
    keep_collinear: bool,
    force: bool,
) -> Option<PolyTree> {
    if !force && (edges.len() < MIN_EDGES || ring_starts.len() < 2) {
        return None;
    }
    let rings = Rings::new(edges, ring_starts);
    let cl = find_clusters(&rings, force)?;
    let nc = cl.start.len() - 1;
    let cluster_rings = |c: usize| &cl.rings[cl.start[c] as usize..cl.start[c + 1] as usize];
    // One query point per cluster: the first vertex of its first ring.
    let kd = KdTree::new(
        (0..nc)
            .map(|c| (rings.edges(cluster_rings(c)[0] as usize)[0].a, c as u32))
            .collect(),
    );
    // Base windings: other clusters' rings around each cluster.
    let mut base = vec![[0i32; 2]; nc];
    let mut hits: Vec<(Point, u32)> = Vec::new();
    let mut budget = 4 * edges.len() + 4096;
    for r in 0..rings.len() {
        let own = cl.of_ring[r];
        hits.clear();
        kd.query(&rings.bbox[r], &mut |q, c| {
            if c != own {
                hits.push((q, c))
            }
        });
        if hits.is_empty() {
            continue;
        }
        budget = match budget.checked_sub(hits.len()) {
            Some(b) => b,
            None if force => budget,
            None => return None,
        };
        let re = rings.edges(r);
        let op = re[0].operand as usize;
        windings(
            re.len(),
            |i| (re[i].a, re[i].b),
            &mut hits,
            |c, w| base[c as usize][op] += w,
        );
    }
    // Every cluster on its own.
    let ids: Vec<u32> = (0..nc as u32).collect();
    let parts: Vec<Part> = crate::par::map_items(&ids, |&c| {
        let rs = cluster_rings(c as usize);
        let b = base[c as usize];
        if rs.len() == 1
            && let Some(p) = lone_ring(rings.edges(rs[0] as usize), b, &inside, keep_collinear)
        {
            return p;
        }
        let mut ce: Vec<InEdge> = Vec::new();
        for &r in rs {
            ce.extend_from_slice(rings.edges(r as usize));
        }
        let dir = compute(&ce, |w| inside([w[0] + b[0], w[1] + b[1]]));
        let (rings, is_hole, parent) = assemble_rings(&dir, keep_collinear);
        Part {
            rings,
            is_hole,
            parent,
        }
    });
    // Parents of the clusters' top-level rings: the innermost output ring of other clusters
    // around them (nested rings never share a vertex, so the smallest area decides).
    let need: Vec<bool> = parts.iter().map(|p| p.parent.contains(&None)).collect();
    let mut offset = Vec::with_capacity(nc + 1);
    offset.push(0u32);
    for p in &parts {
        offset.push(offset.last().unwrap() + p.rings.len() as u32);
    }
    let mut enclosing: Vec<Option<(i128, u32)>> = vec![None; nc];
    for (d, p) in parts.iter().enumerate() {
        for (k, r) in p.rings.iter().enumerate() {
            let Some(bb) = Rect::of_points(&r.pts) else {
                continue;
            };
            hits.clear();
            kd.query(&bb, &mut |q, c| {
                if c as usize != d && need[c as usize] {
                    hits.push((q, c))
                }
            });
            if hits.is_empty() {
                continue;
            }
            let a = ring_area2(&r.pts).abs();
            let g = offset[d] + k as u32;
            let m = r.pts.len();
            windings(
                m,
                |i| (r.pts[i], r.pts[(i + 1) % m]),
                &mut hits,
                |c, _| {
                    let e = &mut enclosing[c as usize];
                    if e.is_none_or(|(a0, _)| a < a0) {
                        *e = Some((a, g));
                    }
                },
            );
        }
    }
    let total = *offset.last().unwrap() as usize;
    let mut all: Vec<RawRing> = Vec::with_capacity(total);
    let mut is_hole: Vec<bool> = Vec::with_capacity(total);
    let mut parent: Vec<Option<u32>> = Vec::with_capacity(total);
    for (c, p) in parts.into_iter().enumerate() {
        let off = offset[c];
        let up = enclosing[c].map(|(_, g)| g);
        parent.extend(p.parent.iter().map(|q| q.map_or(up, |q| Some(q + off))));
        is_hole.extend(p.is_hole);
        all.extend(p.rings);
    }
    Some(canonical_tree(all, &is_hole, &parent))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Quadratic reference for [`ring_is_simple`].
    fn simple_brute(p: &[Point]) -> bool {
        let m = p.len();
        if m < 3 {
            return false;
        }
        for e in 0..m {
            for f in e + 1..m {
                let (a, b) = (p[e], p[(e + 1) % m]);
                let (c, d) = (p[f], p[(f + 1) % m]);
                if !segments_intersect(a, b, c, d) {
                    continue;
                }
                // Adjacent edges: only their common vertex, without folding back.
                let (u, v, w) = if f == e + 1 {
                    (a, b, d)
                } else if e == 0 && f == m - 1 {
                    (c, a, b)
                } else {
                    return false;
                };
                if orient(u, v, w) == 0 && dot(sub(u, v), sub(w, v)) > 0 {
                    return false;
                }
            }
        }
        true
    }

    #[test]
    fn simple_matches_brute_force() {
        let mut s: u64 = 0x5eed;
        let mut rnd = |m: i64| {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((s >> 33) as i64).rem_euclid(m)
        };
        let mut simple = 0;
        for it in 0..20000 {
            let n = 3 + rnd(if it % 3 == 0 { 120 } else { 12 }) as usize;
            let mut pts: Vec<Point> = if it % 2 == 0 {
                // Random points in a small range: crossings, touching, folds.
                let r = 2 + rnd(30);
                (0..n).map(|_| Point::new(rnd(r), rnd(r))).collect()
            } else {
                // Star-shaped (simple unless angles collide).
                let mut a: Vec<(f64, f64)> = (0..n)
                    .map(|_| {
                        let t = rnd(10_000) as f64 / 10_000.0 * core::f64::consts::TAU;
                        (t, 10.0 + rnd(200) as f64)
                    })
                    .collect();
                a.sort_by(|x, y| x.0.total_cmp(&y.0));
                a.iter()
                    .map(|&(t, r)| {
                        Point::new((r * t.cos()).round() as i64, (r * t.sin()).round() as i64)
                    })
                    .collect()
            };
            pts.dedup();
            while pts.len() > 1 && pts[0] == *pts.last().unwrap() {
                pts.pop();
            }
            let (fast, brute) = (ring_is_simple(&pts), simple_brute(&pts));
            assert!(!fast || brute, "{pts:?} reported simple");
            assert!(fast || !brute, "{pts:?} not recognized as simple");
            simple += brute as usize;
        }
        assert!(simple > 500, "{simple}");
    }
}
