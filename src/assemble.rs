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

/// Links boundary edges into rings (as lists of edge indices). Every vertex must have
/// equal in- and out-degree (true for the boundary of any region). At vertices of degree > 2
/// the outgoing edge is chosen first clockwise from the incoming one (face tracing), and
/// closed walks that still visit a vertex twice are split there, so every ring is simple.
/// Also returns the pinch vertices (out-degree > 1), sorted.
pub(crate) fn link_rings(edges: &[DirEdge]) -> (Vec<Vec<u32>>, Vec<Point>) {
    let n = edges.len();
    let mut perm: Vec<u32> = (0..n as u32).collect();
    perm.sort_unstable_by_key(|&k| edges[k as usize].from);
    let from = |pos: usize| edges[perm[pos] as usize].from;
    // Pinch vertices: out-degree > 1.
    let mut pinch: Vec<Point> = Vec::new();
    for pos in 1..n {
        if from(pos) == from(pos - 1) && pinch.last() != Some(&from(pos)) {
            pinch.push(from(pos));
        }
    }
    // First position (in `perm`) of every vertex's out-edges.
    let mut first: PointMap = PointMap::with_capacity(n);
    for pos in 0..n {
        if pos == 0 || from(pos - 1) != from(pos) {
            first.insert(from(pos), pos as u32);
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
            while hi < n && from(hi) == e.to {
                hi += 1;
            }
            let next = if hi - lo == 1 {
                perm[lo] as usize
            } else {
                let r = sub(e.from, e.to);
                let mut best: Option<usize> = None;
                for &k in &perm[lo..hi] {
                    let k = k as usize;
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
            if used[next] {
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
            split_walk(edges, &walk, &mut rings);
        } else {
            rings.push(walk.clone());
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
fn split_walk(edges: &[DirEdge], walk: &[u32], out: &mut Vec<Vec<u32>>) {
    let mut stack: Vec<u32> = Vec::new();
    let mut seen: std::collections::HashMap<Point, usize> = std::collections::HashMap::new();
    for &k in walk {
        let p = edges[k as usize].from;
        if let Some(&pos) = seen.get(&p) {
            let lp: Vec<u32> = stack.drain(pos..).collect();
            for &q in &lp {
                seen.remove(&edges[q as usize].from);
            }
            out.push(lp);
        }
        seen.insert(p, stack.len());
        stack.push(k);
    }
    if !stack.is_empty() {
        out.push(stack);
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

/// Parent of every ring in the nesting tree. `edges` must be in sweep order (as produced
/// by the arrangement), `rings` simple, pairwise non-crossing and oriented with the interior
/// on their left.
fn nesting(edges: &[DirEdge], rings: &[Vec<u32>], is_hole: &[bool]) -> Vec<Option<u32>> {
    let n = edges.len();
    const NONE: u32 = u32::MAX;
    let mut ring_of = vec![NONE; n];
    let mut query = vec![false; n];
    for (ri, r) in rings.iter().enumerate() {
        for &k in r {
            ring_of[k as usize] = ri as u32;
        }
        // At the ring's smallest vertex, the lower of its two edges is the query edge.
        let m = r.len();
        let (i, _) = r
            .iter()
            .enumerate()
            .min_by_key(|(_, k)| edges[**k as usize].from)
            .unwrap();
        let out_e = r[i] as usize;
        let in_e = r[(i + m - 1) % m] as usize;
        let v = edges[out_e].from;
        let d_out = sub(edges[out_e].to, v);
        let d_in = sub(edges[in_e].from, v);
        let q = if cmp_dir_halfplane(d_out, d_in) == Ordering::Less {
            out_e
        } else {
            in_e
        };
        query[q] = true;
    }
    let segs: Vec<(Point, Point)> = edges
        .iter()
        .map(|e| {
            if e.from < e.to {
                (e.from, e.to)
            } else {
                (e.to, e.from)
            }
        })
        .collect();
    debug_assert!(
        segs.windows(2)
            .all(|w| cmp_sweep_edges(w[0], w[1]) != Ordering::Greater)
    );
    let mut parent: Vec<Option<u32>> = vec![None; rings.len()];
    sweep(&segs, |e, below| {
        let e = e as usize;
        if !query[e] || ring_of[e] == NONE {
            return;
        }
        let r = ring_of[e] as usize;
        parent[r] = match below {
            None => None,
            Some(b) => {
                let s = ring_of[b as usize];
                if s == NONE {
                    None
                } else {
                    let s = s as usize;
                    let eb = &edges[b as usize];
                    if eb.from < eb.to {
                        // Material above `b`: r is a hole of the polygon owning s.
                        debug_assert!(is_hole[r]);
                        if is_hole[s] {
                            parent[s]
                        } else {
                            Some(s as u32)
                        }
                    } else {
                        // Exterior (or hole space) above `b`.
                        debug_assert!(!is_hole[r]);
                        if is_hole[s] {
                            Some(s as u32)
                        } else {
                            parent[s]
                        }
                    }
                }
            }
        };
    });
    parent
}

/// Full pipeline from boundary edges (in sweep order) to a canonical tree.
pub(crate) fn assemble(edges: Vec<DirEdge>, keep_collinear: bool) -> PolyTree {
    let (rings_idx, pinch) = link_rings(&edges);
    let pts_of =
        |r: &Vec<u32>| -> Vec<Point> { r.iter().map(|&k| edges[k as usize].from).collect() };
    // Degenerate loops cannot occur for valid arrangements; drop them defensively.
    let rings_idx: Vec<Vec<u32>> = rings_idx
        .into_iter()
        .filter(|r| r.len() >= 3 && ring_area2(&pts_of(r)) != 0)
        .collect();
    let is_hole: Vec<bool> = rings_idx
        .iter()
        .map(|r| ring_area2(&pts_of(r)) < 0)
        .collect();
    let parent = nesting(&edges, &rings_idx, &is_hole);
    let mut rings: Vec<RawRing> = rings_idx
        .iter()
        .map(|r| RawRing {
            pts: pts_of(r),
            tags: r.iter().map(|&k| edges[k as usize].tag).collect(),
        })
        .collect();
    drop(rings_idx);
    for r in rings.iter_mut() {
        if !keep_collinear {
            remove_collinear(r, &pinch);
        }
        rotate_to_min(r);
    }
    canonical_tree(rings, &is_hole, &parent)
}

/// Orders the rings canonically (children by vertex sequence, depth-first numbering).
fn canonical_tree(rings: Vec<RawRing>, is_hole: &[bool], parent: &[Option<u32>]) -> PolyTree {
    let m = rings.len();
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
