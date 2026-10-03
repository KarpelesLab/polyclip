//! Noding: splitting a set of segments at all their mutual intersections.
//!
//! Two modes are provided:
//!
//! * [`snap_round`] — Hobby's snap rounding. Every segment endpoint and every rounded proper
//!   crossing defines a *hot pixel* (the half-open unit square centred on an integer point).
//!   Each segment is replaced by the polyline through the centres of the hot pixels it meets,
//!   in the order it meets them. The result has integer vertices, no proper crossings, and
//!   moves no point by more than `sqrt(2)/2`. A final pass splits fragments at any hot pixel
//!   centre lying exactly on their interior, so the output is fully noded: fragments meet only
//!   at shared endpoints or coincide exactly.
//! * [`node_exact`] — no rounding at all: fails if any two segments cross properly, otherwise
//!   splits segments at vertices lying on their interiors.
//!
//! Both use a k-d partition of the plane into leaves holding few segments each.

use crate::geom::{Point, Rect};
use crate::predicates::{
    Frac, in_segment_interior, orient, rounded_crossing, segment_meets_rect, segment_pixel_entry,
    segments_cross_properly,
};

/// A fragment of input segment `src`, from `a` to `b` (same direction as the source).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Frag {
    pub a: Point,
    pub b: Point,
    pub src: u32,
}

const LEAF_SIZE: usize = 16;
const MAX_DEPTH: u32 = 48;

#[derive(Clone, Copy, Debug)]
enum KdNode {
    Split {
        axis: u8,
        at: i64,
        left: u32,
        right: u32,
    },
    Leaf {
        leaf: u32,
    },
}

/// A leaf: a half-open integer region `[lo, hi)` and its segments.
#[derive(Clone, Debug)]
struct Leaf {
    lo: Point,
    hi: Point,
    start: u32,
    end: u32,
}

struct Kd {
    nodes: Vec<KdNode>,
    leaves: Vec<Leaf>,
    items: Vec<u32>,
}

#[inline]
fn seg_bbox(s: &(Point, Point)) -> Rect {
    Rect::new(s.0, s.1)
}

impl Kd {
    /// Builds the partition. A segment belongs to every leaf whose region, grown by one unit,
    /// it meets; this guarantees that a segment passing through (or within distance 1 of) an
    /// integer point of a leaf region is listed in that leaf.
    fn build(segs: &[(Point, Point)], bboxes: &[Rect]) -> Kd {
        let mut kd = Kd {
            nodes: Vec::new(),
            leaves: Vec::new(),
            items: Vec::new(),
        };
        let Some(root_box) = bboxes.iter().copied().reduce(|a, b| a.union(&b)) else {
            return kd;
        };
        let lo = root_box.min;
        let hi = Point::new(root_box.max.x + 1, root_box.max.y + 1);
        kd.nodes.push(KdNode::Leaf { leaf: 0 });
        let all: Vec<u32> = (0..segs.len() as u32).collect();
        let mut stack: Vec<(u32, Vec<u32>, Point, Point, u32)> = vec![(0, all, lo, hi, 0)];
        let mut centers: Vec<i64> = Vec::new();
        while let Some((node, items, lo, hi, depth)) = stack.pop() {
            let split = if items.len() > LEAF_SIZE && depth < MAX_DEPTH {
                Self::choose_split(&items, bboxes, lo, hi, &mut centers)
            } else {
                None
            };
            if let Some((axis, at)) = split {
                let (lhi, rlo) = if axis == 0 {
                    (Point::new(at, hi.y), Point::new(at, lo.y))
                } else {
                    (Point::new(hi.x, at), Point::new(lo.x, at))
                };
                let lrect = Rect {
                    min: Point::new(lo.x - 1, lo.y - 1),
                    max: lhi,
                };
                let rrect = Rect {
                    min: Point::new(rlo.x - 1, rlo.y - 1),
                    max: hi,
                };
                let mut left = Vec::new();
                let mut right = Vec::new();
                for &i in &items {
                    let s = &segs[i as usize];
                    let b = &bboxes[i as usize];
                    if b.intersects(&lrect) && segment_meets_rect(s.0, s.1, &lrect) {
                        left.push(i);
                    }
                    if b.intersects(&rrect) && segment_meets_rect(s.0, s.1, &rrect) {
                        right.push(i);
                    }
                }
                // Stop when splitting mostly duplicates items (long segments).
                if (left.len() + right.len()) * 10 <= items.len() * 19 {
                    let l = kd.nodes.len() as u32;
                    kd.nodes.push(KdNode::Leaf { leaf: 0 });
                    kd.nodes.push(KdNode::Leaf { leaf: 0 });
                    kd.nodes[node as usize] = KdNode::Split {
                        axis,
                        at,
                        left: l,
                        right: l + 1,
                    };
                    stack.push((l + 1, right, rlo, hi, depth + 1));
                    stack.push((l, left, lo, lhi, depth + 1));
                    continue;
                }
            }
            let leaf = kd.leaves.len() as u32;
            let start = kd.items.len() as u32;
            kd.items.extend_from_slice(&items);
            kd.leaves.push(Leaf {
                lo,
                hi,
                start,
                end: kd.items.len() as u32,
            });
            kd.nodes[node as usize] = KdNode::Leaf { leaf };
        }
        kd
    }

    fn choose_split(
        items: &[u32],
        bboxes: &[Rect],
        lo: Point,
        hi: Point,
        centers: &mut Vec<i64>,
    ) -> Option<(u8, i64)> {
        let w = hi.x - lo.x;
        let h = hi.y - lo.y;
        let axes: [u8; 2] = if w >= h { [0, 1] } else { [1, 0] };
        for axis in axes {
            let (l, r) = if axis == 0 {
                (lo.x, hi.x)
            } else {
                (lo.y, hi.y)
            };
            if r - l < 2 {
                continue;
            }
            centers.clear();
            centers.extend(items.iter().map(|&i| {
                let b = &bboxes[i as usize];
                if axis == 0 {
                    b.min.x + (b.max.x - b.min.x) / 2
                } else {
                    b.min.y + (b.max.y - b.min.y) / 2
                }
            }));
            let mid = centers.len() / 2;
            let (_, m, _) = centers.select_nth_unstable(mid);
            let mut at = *m;
            if at <= l || at >= r {
                at = l + (r - l) / 2;
            }
            return Some((axis, at.clamp(l + 1, r - 1)));
        }
        None
    }

    /// Index of the leaf whose region contains `p` (which must lie in the root region).
    fn locate(&self, p: Point) -> u32 {
        let mut n = 0usize;
        loop {
            match self.nodes[n] {
                KdNode::Leaf { leaf } => return leaf,
                KdNode::Split {
                    axis,
                    at,
                    left,
                    right,
                } => {
                    let c = if axis == 0 { p.x } else { p.y };
                    n = if c < at {
                        left as usize
                    } else {
                        right as usize
                    };
                }
            }
        }
    }

    fn leaf_items(&self, l: &Leaf) -> &[u32] {
        &self.items[l.start as usize..l.end as usize]
    }
}

#[inline]
fn in_region(l: &Leaf, p: Point) -> bool {
    p.x >= l.lo.x && p.x < l.hi.x && p.y >= l.lo.y && p.y < l.hi.y
}

/// Calls `f(i, j)` for every pair of segments in the leaf whose bounding boxes overlap.
fn for_each_pair(
    items: &[u32],
    bboxes: &[Rect],
    order: &mut Vec<u32>,
    mut f: impl FnMut(u32, u32),
) {
    order.clear();
    order.extend_from_slice(items);
    order.sort_unstable_by_key(|&i| bboxes[i as usize].min.x);
    for (k, &i) in order.iter().enumerate() {
        let bi = &bboxes[i as usize];
        for &j in &order[k + 1..] {
            let bj = &bboxes[j as usize];
            if bj.min.x > bi.max.x {
                break;
            }
            if bj.min.y <= bi.max.y && bi.min.y <= bj.max.y {
                f(i, j);
            }
        }
    }
}

/// Snap-rounds `segs` (each with distinct endpoints). Returns the fragments of every
/// segment, grouped by segment in input order and ordered from `a` to `b` within a segment.
pub(crate) fn snap_round(segs: &[(Point, Point)]) -> Vec<Frag> {
    let bboxes: Vec<Rect> = segs.iter().map(seg_bbox).collect();
    let kd = Kd::build(segs, &bboxes);
    let mut order = Vec::new();

    // 1. Hot pixels: endpoints and rounded proper crossings.
    let mut hot: Vec<Point> = Vec::with_capacity(segs.len() * 2);
    for s in segs {
        hot.push(s.0);
        hot.push(s.1);
    }
    for leaf in &kd.leaves {
        for_each_pair(kd.leaf_items(leaf), &bboxes, &mut order, |i, j| {
            let (a, b) = segs[i as usize];
            let (c, d) = segs[j as usize];
            if segments_cross_properly(a, b, c, d) {
                let p = rounded_crossing(a, b, c, d);
                // Report each crossing once: in the leaf containing its pixel.
                if in_region(leaf, p) {
                    hot.push(p);
                }
            }
        });
    }
    hot.sort_unstable();
    hot.dedup();

    // 2. Distribute hot pixels to leaves, sorted by (leaf, x, y).
    let mut by_leaf: Vec<(u32, Point)> = hot.iter().map(|&p| (kd.locate(p), p)).collect();
    drop(hot);
    by_leaf.sort_unstable();
    let mut leaf_start = vec![0u32; kd.leaves.len() + 1];
    for &(l, _) in &by_leaf {
        leaf_start[l as usize + 1] += 1;
    }
    for i in 0..kd.leaves.len() {
        leaf_start[i + 1] += leaf_start[i];
    }

    // 3. For each segment, the hot pixels it meets (other than its own endpoints) and the
    //    hot pixel centres close to it.
    struct Hit {
        seg: u32,
        t: Frac,
        open: bool,
        p: Point,
    }
    let mut hits: Vec<Hit> = Vec::new();
    let mut near: Vec<(u32, Point)> = Vec::new();
    for (li, leaf) in kd.leaves.iter().enumerate() {
        let pix = &by_leaf[leaf_start[li] as usize..leaf_start[li + 1] as usize];
        if pix.is_empty() {
            continue;
        }
        for &si in kd.leaf_items(leaf) {
            let (a, b) = segs[si as usize];
            let bb = &bboxes[si as usize];
            let x0 = bb.min.x - 1;
            let x1 = bb.max.x + 1;
            let from = pix.partition_point(|&(_, p)| p.x < x0);
            let len2 = {
                let dx = (b.x - a.x) as f64;
                let dy = (b.y - a.y) as f64;
                dx * dx + dy * dy
            };
            for &(_, p) in &pix[from..] {
                if p.x > x1 {
                    break;
                }
                if p.y < bb.min.y - 1 || p.y > bb.max.y + 1 || p == a || p == b {
                    continue;
                }
                if let Some((t, open)) = segment_pixel_entry(a, b, p) {
                    hits.push(Hit {
                        seg: si,
                        t,
                        open,
                        p,
                    });
                } else {
                    // Distance filter (approximate, conservative): only centres within ~1
                    // unit of the segment can lie on one of its fragments.
                    let o = orient(a, b, p) as f64;
                    if o * o <= 4.0 * len2 {
                        near.push((si, p));
                    }
                }
            }
        }
    }
    drop(by_leaf);
    // Pixels partition the plane, so (entry value, closed-before-open) is a strict order.
    hits.sort_unstable_by(|x, y| {
        x.seg
            .cmp(&y.seg)
            .then_with(|| x.t.cmp(y.t))
            .then(x.open.cmp(&y.open))
    });
    near.sort_unstable();

    // 4. Build fragments.
    let mut out: Vec<Frag> = Vec::with_capacity(segs.len() + hits.len());
    let mut poly: Vec<Point> = Vec::new();
    let mut ins: Vec<(usize, i128, Point)> = Vec::new();
    let mut hi = 0usize;
    let mut ni = 0usize;
    for (si, &(a, b)) in segs.iter().enumerate() {
        let si = si as u32;
        poly.clear();
        poly.push(a);
        while hi < hits.len() && hits[hi].seg == si {
            let p = hits[hi].p;
            if *poly.last().unwrap() != p {
                poly.push(p);
            }
            hi += 1;
        }
        if *poly.last().unwrap() != b {
            poly.push(b);
        }
        // T-junctions: hot pixel centres lying on a fragment's interior.
        ins.clear();
        while ni < near.len() && near[ni].0 == si {
            let c = near[ni].1;
            ni += 1;
            for k in 0..poly.len() - 1 {
                if in_segment_interior(poly[k], poly[k + 1], c) {
                    ins.push((k, crate::predicates::dist2(poly[k], c), c));
                    break;
                }
            }
        }
        if ins.is_empty() {
            for w in poly.windows(2) {
                out.push(Frag {
                    a: w[0],
                    b: w[1],
                    src: si,
                });
            }
        } else {
            ins.sort_unstable();
            let mut j = 0;
            for k in 0..poly.len() - 1 {
                let mut cur = poly[k];
                while j < ins.len() && ins[j].0 == k {
                    if ins[j].2 != cur {
                        out.push(Frag {
                            a: cur,
                            b: ins[j].2,
                            src: si,
                        });
                        cur = ins[j].2;
                    }
                    j += 1;
                }
                out.push(Frag {
                    a: cur,
                    b: poly[k + 1],
                    src: si,
                });
            }
        }
    }
    out
}

/// A proper crossing found by [`node_exact`]: the two segment indices (`i < j`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct Crossing {
    pub i: u32,
    pub j: u32,
}

/// Exact noding without rounding. Returns the smallest crossing pair if any two segments
/// cross properly; otherwise splits each segment at every segment endpoint lying on its
/// interior and returns the fragments (grouped by segment, ordered from `a` to `b`).
pub(crate) fn node_exact(segs: &[(Point, Point)]) -> Result<Vec<Frag>, Crossing> {
    let bboxes: Vec<Rect> = segs.iter().map(seg_bbox).collect();
    let kd = Kd::build(segs, &bboxes);
    let mut order = Vec::new();
    let mut crossing: Option<Crossing> = None;
    let mut splits: Vec<(u32, i128, Point)> = Vec::new();
    for leaf in &kd.leaves {
        for_each_pair(kd.leaf_items(leaf), &bboxes, &mut order, |i, j| {
            let (a, b) = segs[i as usize];
            let (c, d) = segs[j as usize];
            if segments_cross_properly(a, b, c, d) {
                let cr = Crossing {
                    i: i.min(j),
                    j: i.max(j),
                };
                if crossing.is_none_or(|c0| cr < c0) {
                    crossing = Some(cr);
                }
                return;
            }
            for (s, (p, q), pts) in [(i, (a, b), [c, d]), (j, (c, d), [a, b])] {
                for v in pts {
                    if in_region(leaf, v) && in_segment_interior(p, q, v) {
                        splits.push((s, crate::predicates::dist2(p, v), v));
                    }
                }
            }
        });
    }
    if let Some(c) = crossing {
        return Err(c);
    }
    splits.sort_unstable();
    splits.dedup();
    let mut out = Vec::with_capacity(segs.len() + splits.len());
    let mut k = 0;
    for (si, &(a, b)) in segs.iter().enumerate() {
        let si = si as u32;
        let mut cur = a;
        while k < splits.len() && splits[k].0 == si {
            out.push(Frag {
                a: cur,
                b: splits[k].2,
                src: si,
            });
            cur = splits[k].2;
            k += 1;
        }
        out.push(Frag { a: cur, b, src: si });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::predicates::segments_cross_properly;

    fn p(x: i64, y: i64) -> Point {
        Point::new(x, y)
    }

    /// Checks that fragments are fully noded: no proper crossings, no vertex on another
    /// fragment's interior.
    pub(crate) fn assert_noded(frags: &[Frag]) {
        for (i, f) in frags.iter().enumerate() {
            assert_ne!(f.a, f.b);
            for g in &frags[i + 1..] {
                assert!(
                    !segments_cross_properly(f.a, f.b, g.a, g.b),
                    "{f:?} crosses {g:?}"
                );
                for v in [g.a, g.b] {
                    assert!(!in_segment_interior(f.a, f.b, v), "{v:?} on {f:?}");
                }
                for v in [f.a, f.b] {
                    assert!(!in_segment_interior(g.a, g.b, v), "{v:?} on {g:?}");
                }
            }
        }
    }

    #[test]
    fn simple_cross() {
        let segs = [(p(0, 0), p(10, 10)), (p(0, 10), p(10, 0))];
        let f = snap_round(&segs);
        assert_eq!(f.len(), 4);
        assert!(f.iter().all(|f| f.a == p(5, 5) || f.b == p(5, 5)));
        assert_noded(&f);
    }

    #[test]
    fn rounded_cross() {
        let segs = [(p(0, 0), p(3, 3)), (p(0, 3), p(3, 0))];
        let f = snap_round(&segs);
        // Crossing at (1.5, 1.5) -> pixel (2, 2).
        assert!(f.iter().any(|f| f.b == p(2, 2)));
        assert_noded(&f);
    }

    #[test]
    fn exact_t_junction() {
        let segs = [(p(0, 0), p(10, 0)), (p(5, 0), p(5, 5))];
        let f = node_exact(&segs).unwrap();
        assert_eq!(f.len(), 3);
        assert_noded(&f);
        let segs = [(p(0, 0), p(10, 0)), (p(5, -1), p(5, 5))];
        assert!(node_exact(&segs).is_err());
    }

    #[test]
    fn random_snap_noded() {
        // Deterministic LCG for a quick self-contained stress test.
        let mut s: u64 = 0x1234_5678;
        let mut rnd = |m: i64| {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((s >> 33) as i64).rem_euclid(m)
        };
        for it in 0..600 {
            let n = 2 + rnd(60) as usize;
            let range = 1 + rnd(if it % 2 == 0 { 30 } else { 3000 });
            let segs: Vec<(Point, Point)> = (0..n)
                .map(|_| (p(rnd(range), rnd(range)), p(rnd(range), rnd(range))))
                .filter(|(a, b)| a != b)
                .collect();
            let f = snap_round(&segs);
            assert_noded(&f);
        }
    }
}
