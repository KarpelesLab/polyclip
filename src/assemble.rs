//! Ring assembly: turns the directed boundary edges of a result region (interior on the
//! left of every edge) into canonical rings and their nesting tree.

use crate::geom::{Point, PolyNode, PolyTree, Ring};
use crate::predicates::{cmp_dir_halfplane, cross, dot, orient, sub};
use crate::query::ring_area2;
use crate::sweep::{cmp_sweep_edges, sweep};
use core::cmp::Ordering;

/// A directed boundary edge with its provenance tag.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct DirEdge {
    pub from: Point,
    pub to: Point,
    pub tag: u64,
}

/// A ring under construction: vertices and per-edge tags.
#[derive(Clone, Debug)]
pub(crate) struct RawRing {
    pub pts: Vec<Point>,
    pub tags: Vec<u64>,
}

/// Compares the clockwise angle from `r` to `a` with the clockwise angle from `r` to `b`
/// (all non-zero vectors, neither `a` nor `b` pointing along `r`).
fn cmp_cw_from(r: Point, a: Point, b: Point) -> Ordering {
    // Half 0: clockwise angle in (0, pi]; half 1: (pi, 2pi).
    let half = |v: Point| {
        let c = cross(r, v);
        if c < 0 || (c == 0 && dot(r, v) < 0) {
            0u8
        } else {
            1u8
        }
    };
    half(a).cmp(&half(b)).then_with(|| cross(a, b).cmp(&0))
}

/// Links boundary edges into rings. Every vertex must have equal in- and out-degree (true
/// for the boundary of any region). At vertices of degree > 2 the outgoing edge is chosen
/// first clockwise from the incoming one (face tracing), and closed walks that still visit a
/// vertex twice are split there, so every returned ring is simple.
pub(crate) fn link_rings(mut edges: Vec<DirEdge>) -> (Vec<RawRing>, Vec<Point>) {
    edges.sort_unstable();
    let n = edges.len();
    // Pinch vertices: out-degree > 1.
    let mut pinch: Vec<Point> = Vec::new();
    for w in edges.windows(2) {
        if w[0].from == w[1].from && pinch.last() != Some(&w[0].from) {
            pinch.push(w[0].from);
        }
    }
    // First out-edge index of every vertex.
    let mut first: PointMap = PointMap::with_capacity(n);
    for (i, e) in edges.iter().enumerate() {
        if i == 0 || edges[i - 1].from != e.from {
            first.insert(e.from, i as u32);
        }
    }
    let mut used = vec![false; n];
    let mut rings = Vec::new();
    let mut walk: Vec<u32> = Vec::new();
    for start in 0..n {
        if used[start] {
            continue;
        }
        walk.clear();
        let mut cur = start;
        let mut closed = false;
        loop {
            used[cur] = true;
            walk.push(cur as u32);
            let e = edges[cur];
            let Some(lo) = first.get(e.to) else { break };
            let lo = lo as usize;
            let mut hi = lo + 1;
            while hi < n && edges[hi].from == e.to {
                hi += 1;
            }
            let next = if hi - lo == 1 {
                lo
            } else {
                let r = sub(e.from, e.to);
                let mut best: Option<usize> = None;
                for k in lo..hi {
                    if used[k] && k != start {
                        continue;
                    }
                    let d = sub(edges[k].to, edges[k].from);
                    if best.is_none_or(|b| {
                        cmp_cw_from(r, d, sub(edges[b].to, edges[b].from)) == Ordering::Less
                    }) {
                        best = Some(k);
                    }
                }
                match best {
                    Some(b) => b,
                    None => break,
                }
            };
            if next == start {
                closed = true;
                break;
            }
            if hi == lo || used[next] {
                break;
            }
            cur = next;
        }
        if !closed {
            debug_assert!(false, "link_rings: open walk");
            continue;
        }
        let has_pinch = !pinch.is_empty()
            && walk
                .iter()
                .any(|&k| pinch.binary_search(&edges[k as usize].from).is_ok());
        if has_pinch {
            split_walk(&edges, &walk, &mut rings);
        } else {
            rings.push(RawRing {
                pts: walk.iter().map(|&k| edges[k as usize].from).collect(),
                tags: walk.iter().map(|&k| edges[k as usize].tag).collect(),
            });
        }
    }
    (rings, pinch)
}

/// A minimal open-addressing hash map from points to `u32`, with a fixed multiplicative
/// hash (lookups only; iteration order never matters, so results stay deterministic).
pub(crate) struct PointMap {
    keys: Vec<Point>,
    vals: Vec<u32>,
    mask: usize,
}

impl PointMap {
    const EMPTY: u32 = u32::MAX;

    pub fn with_capacity(n: usize) -> Self {
        let cap = (n * 2).next_power_of_two().max(16);
        PointMap {
            keys: vec![Point::default(); cap],
            vals: vec![Self::EMPTY; cap],
            mask: cap - 1,
        }
    }

    #[inline]
    fn slot(&self, p: Point) -> usize {
        let h = (p.x as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15)
            ^ (p.y as u64).wrapping_mul(0xC2B2_AE3D_27D4_EB4F);
        ((h ^ (h >> 29)).wrapping_mul(0x1656_67B1_9E37_79F9) >> 20) as usize & self.mask
    }

    /// Inserts (or overwrites) `p -> v`; `v` must not be `u32::MAX`.
    pub fn insert(&mut self, p: Point, v: u32) {
        let mut i = self.slot(p);
        loop {
            if self.vals[i] == Self::EMPTY || self.keys[i] == p {
                self.keys[i] = p;
                self.vals[i] = v;
                return;
            }
            i = (i + 1) & self.mask;
        }
    }

    #[inline]
    pub fn get(&self, p: Point) -> Option<u32> {
        let mut i = self.slot(p);
        loop {
            let v = self.vals[i];
            if v == Self::EMPTY {
                return None;
            }
            if self.keys[i] == p {
                return Some(v);
            }
            i = (i + 1) & self.mask;
        }
    }
}

/// Splits a closed walk that may repeat vertices into simple loops.
fn split_walk(edges: &[DirEdge], walk: &[u32], out: &mut Vec<RawRing>) {
    let mut stack: Vec<u32> = Vec::new();
    let mut seen: std::collections::HashMap<Point, usize> = std::collections::HashMap::new();
    for &k in walk {
        let p = edges[k as usize].from;
        if let Some(&pos) = seen.get(&p) {
            let lp: Vec<u32> = stack.drain(pos..).collect();
            for &q in &lp {
                seen.remove(&edges[q as usize].from);
            }
            out.push(RawRing {
                pts: lp.iter().map(|&q| edges[q as usize].from).collect(),
                tags: lp.iter().map(|&q| edges[q as usize].tag).collect(),
            });
        }
        seen.insert(p, stack.len());
        stack.push(k);
    }
    if !stack.is_empty() {
        out.push(RawRing {
            pts: stack.iter().map(|&q| edges[q as usize].from).collect(),
            tags: stack.iter().map(|&q| edges[q as usize].tag).collect(),
        });
    }
}

/// Removes vertices that lie on the straight line through their neighbours, unless they are
/// pinch vertices (shared with other rings) or separate edges with different tags.
pub(crate) fn remove_collinear(r: &mut RawRing, pinch: &[Point]) {
    let n = r.pts.len();
    if n < 3 {
        return;
    }
    let removable = |prev: Point, v: Point, next: Point, tin: u64, tout: u64| {
        tin == tout
            && orient(prev, v, next) == 0
            && (pinch.is_empty() || pinch.binary_search(&v).is_err())
    };
    // Start from a vertex that is kept, so the wrap-around needs no special handling.
    let Some(s) = (0..n).find(|&i| {
        !removable(
            r.pts[(i + n - 1) % n],
            r.pts[i],
            r.pts[(i + 1) % n],
            r.tags[(i + n - 1) % n],
            r.tags[i],
        )
    }) else {
        return;
    };
    let mut pts = Vec::with_capacity(n);
    let mut tags = Vec::with_capacity(n);
    pts.push(r.pts[s]);
    tags.push(r.tags[s]);
    for k in 1..n {
        let i = (s + k) % n;
        let next = r.pts[(i + 1) % n];
        let prev = *pts.last().unwrap();
        if removable(prev, r.pts[i], next, *tags.last().unwrap(), r.tags[i]) {
            continue;
        }
        pts.push(r.pts[i]);
        tags.push(r.tags[i]);
    }
    r.pts = pts;
    r.tags = tags;
}

/// Rotates the ring to start at its lexicographically smallest vertex.
pub(crate) fn rotate_to_min(r: &mut RawRing) {
    if let Some((i, _)) = r.pts.iter().enumerate().min_by_key(|(_, p)| **p) {
        r.pts.rotate_left(i);
        r.tags.rotate_left(i);
    }
}

/// Builds the canonical nesting tree from simple, pairwise non-crossing rings that are
/// oriented with the interior on their left (outer rings counter-clockwise, holes clockwise).
pub(crate) fn build_tree(mut rings: Vec<RawRing>) -> PolyTree {
    for r in rings.iter_mut() {
        rotate_to_min(r);
    }
    let m = rings.len();
    let is_hole: Vec<bool> = rings.iter().map(|r| ring_area2(&r.pts) < 0).collect();

    // Sweep edges: (lo, hi), ring id, material above, query flag.
    struct E {
        lo: Point,
        hi: Point,
        ring: u32,
        inside_above: bool,
        query: bool,
    }
    let mut es: Vec<E> = Vec::with_capacity(rings.iter().map(|r| r.pts.len()).sum());
    for (ri, r) in rings.iter().enumerate() {
        let n = r.pts.len();
        let v = r.pts[0];
        // The lower of the two edges at the minimum vertex is the query edge.
        let d_next = sub(r.pts[1 % n], v);
        let d_prev = sub(r.pts[n - 1], v);
        let next_is_lower = cmp_dir_halfplane(d_next, d_prev) == Ordering::Less;
        for i in 0..n {
            let a = r.pts[i];
            let b = r.pts[(i + 1) % n];
            let (lo, hi) = if a < b { (a, b) } else { (b, a) };
            let query = (i == 0 && next_is_lower) || (i == n - 1 && !next_is_lower);
            es.push(E {
                lo,
                hi,
                ring: ri as u32,
                inside_above: a < b,
                query,
            });
        }
    }
    es.sort_unstable_by(|a, b| cmp_sweep_edges((a.lo, a.hi), (b.lo, b.hi)));
    let segs: Vec<(Point, Point)> = es.iter().map(|e| (e.lo, e.hi)).collect();
    let mut parent: Vec<Option<u32>> = vec![None; m];
    sweep(&segs, |e, below| {
        let e = &es[e as usize];
        if !e.query {
            return;
        }
        let r = e.ring as usize;
        parent[r] = match below {
            None => None,
            Some(b) => {
                let b = &es[b as usize];
                let s = b.ring as usize;
                if b.inside_above {
                    // Material between: r is a hole of the polygon owning s.
                    debug_assert!(is_hole[r]);
                    if is_hole[s] {
                        parent[s]
                    } else {
                        Some(s as u32)
                    }
                } else {
                    // Exterior (or hole space) between.
                    debug_assert!(!is_hole[r]);
                    if is_hole[s] {
                        Some(s as u32)
                    } else {
                        parent[s]
                    }
                }
            }
        };
    });

    // Canonical order: children sorted by ring vertex sequence; nodes numbered depth-first.
    let mut children: Vec<Vec<u32>> = vec![Vec::new(); m];
    let mut roots: Vec<u32> = Vec::new();
    for (r, par) in parent.iter().enumerate() {
        match par {
            Some(p) => children[*p as usize].push(r as u32),
            None => roots.push(r as u32),
        }
    }
    let key = |a: &u32, b: &u32| rings[*a as usize].pts.cmp(&rings[*b as usize].pts);
    roots.sort_by(key);
    for c in children.iter_mut() {
        c.sort_by(key);
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
    PolyTree {
        nodes,
        roots: roots.iter().map(|&r| new_id[r as usize]).collect(),
    }
}

/// Full pipeline from boundary edges to a canonical tree.
pub(crate) fn assemble(edges: Vec<DirEdge>, keep_collinear: bool) -> PolyTree {
    let (mut rings, pinch) = link_rings(edges);
    if !keep_collinear {
        for r in rings.iter_mut() {
            remove_collinear(r, &pinch);
        }
    }
    rings.retain(|r| r.pts.len() >= 3 && ring_area2(&r.pts) != 0);
    build_tree(rings)
}
