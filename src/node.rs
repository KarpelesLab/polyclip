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
//! Both use a uniform grid whose cell size follows the typical segment length; a segment is
//! listed in every cell it passes within distance 1 of, so any two segments meeting near a
//! point are both listed in that point's cell.

use crate::dir::Dir;
use crate::geom::{Point, Rect};
use crate::predicates::{
    dot, floor_div, in_segment_interior, orient, segment_meets_rect, segment_pixel_entry,
    segments_cross_properly, sub,
};

/// A fragment of input segment `src`, from `a` to `b` (same direction as the source).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Frag {
    pub a: Point,
    pub b: Point,
    pub src: u32,
}

/// Adaptive k-d partition of the plane into leaves with half-open integer regions
/// `[lo, hi)`. A segment is listed in every leaf whose region, grown by one unit on every
/// side, it meets — so any two segments passing within distance 1 of an integer point of a
/// leaf are both listed in that leaf. Splits adapt to the data (long, dense or parallel
/// segments included), which a uniform grid cannot do.
struct Grid {
    leaves: Vec<(Point, Point)>,
    start: Vec<u32>,
    items: Vec<u32>,
}

const LEAF_SIZE: usize = 48;
const MAX_DEPTH: u32 = 48;

#[inline]
fn seg_bbox(s: &(Point, Point)) -> Rect {
    Rect::new(s.0, s.1)
}

/// Counting sort of `(key, value)` pairs by key into CSR form. Stable.
fn csr<T: Copy + Default>(n_keys: usize, pairs: &[(u32, T)]) -> (Vec<u32>, Vec<T>) {
    let mut start = vec![0u32; n_keys + 1];
    for &(k, _) in pairs {
        start[k as usize + 1] += 1;
    }
    for i in 0..n_keys {
        start[i + 1] += start[i];
    }
    let mut pos: Vec<u32> = start[..n_keys].to_vec();
    let mut out = vec![T::default(); pairs.len()];
    for &(k, v) in pairs {
        let p = &mut pos[k as usize];
        out[*p as usize] = v;
        *p += 1;
    }
    (start, out)
}

#[inline]
fn coord(p: Point, axis: u8) -> i64 {
    if axis == 0 { p.x } else { p.y }
}

impl Grid {
    /// A uniform grid when the data suits one (segments short relative to their spacing,
    /// little duplication), else an adaptive k-d partition.
    fn build(segs: &[(Point, Point)], bboxes: &[Rect]) -> Grid {
        Self::build_uniform(segs, bboxes).unwrap_or_else(|| Self::build_kd(segs, bboxes))
    }

    /// Uniform grid of square cells sized after the typical segment extent. Returns `None`
    /// when segments would be listed in too many cells (long or dense parallel segments).
    fn build_uniform(segs: &[(Point, Point)], bboxes: &[Rect]) -> Option<Grid> {
        let bb = bboxes.iter().copied().reduce(|a, b| a.union(&b))?;
        let n = segs.len();
        let w = (bb.max.x - bb.min.x + 1) as f64;
        let h = (bb.max.y - bb.min.y + 1) as f64;
        // Typical segment extent (median of a deterministic sample).
        let stride = (n / 4096).max(1);
        let mut ext: Vec<i64> = bboxes
            .iter()
            .step_by(stride)
            .map(|b| b.width().max(b.height()))
            .collect();
        let mid = ext.len() / 2;
        let typical = *ext.select_nth_unstable(mid).1;
        // At most about one cell per segment.
        let s_min = libm::ceil(libm::sqrt(w * h / (n as f64 + 16.0)));
        let mut s = (typical as f64).max(s_min).max(1.0).min(w.max(h)) as i64;
        // Bound the number of cells too (a thin, wide box would otherwise need a huge
        // grid): at most about four cells per segment.
        let cells = |s: i64| {
            ((bb.max.x - bb.min.x) / s + 1) as u128 * ((bb.max.y - bb.min.y) / s + 1) as u128
        };
        while cells(s) > 4 * n as u128 + 64 {
            s = s.saturating_mul(2);
        }
        let nx = ((bb.max.x - bb.min.x) / s + 1) as usize;
        let ny = ((bb.max.y - bb.min.y) / s + 1) as usize;
        let (x0, y0) = (bb.min.x, bb.min.y);
        let col = |x: i64| ((x - x0).div_euclid(s)).clamp(0, nx as i64 - 1) as usize;
        let row = |y: i64| ((y - y0).div_euclid(s)).clamp(0, ny as i64 - 1) as usize;
        // Cells of segment `i`: (columns, rows per column) with the long-diagonal refinement.
        let cells = |i: usize, f: &mut dyn FnMut(usize, usize)| {
            let (sg, b) = (&segs[i], &bboxes[i]);
            let (cx0, cx1) = (col(b.min.x - 1), col(b.max.x + 1));
            let (cy0, cy1) = (row(b.min.y - 1), row(b.max.y + 1));
            // Short or axis-parallel segments: the whole (grown) bounding box of cells.
            if cx1 - cx0 <= 2 || cy1 - cy0 <= 2 || sg.0.x == sg.1.x || sg.0.y == sg.1.y {
                for cy in cy0..=cy1 {
                    for cx in cx0..=cx1 {
                        f(cx, cy);
                    }
                }
                return;
            }
            // Long diagonal segment: per column, the rows its (grown) trace covers.
            let (a, c) = if sg.0.x <= sg.1.x {
                (sg.0, sg.1)
            } else {
                (sg.1, sg.0)
            };
            let slope = (c.y - a.y) as f64 / (c.x - a.x) as f64;
            for cx in cx0..=cx1 {
                let xa = (x0 + cx as i64 * s - 1).clamp(a.x, c.x);
                let xb = (x0 + (cx as i64 + 1) * s).clamp(a.x, c.x);
                let ya = a.y as f64 + (xa - a.x) as f64 * slope;
                let yb = a.y as f64 + (xb - a.x) as f64 * slope;
                let lo = libm::floor(ya.min(yb)) as i64 - 2;
                let hi = libm::ceil(ya.max(yb)) as i64 + 2;
                for cy in row(lo).max(cy0)..=row(hi).min(cy1) {
                    f(cx, cy);
                }
            }
        };
        // Budget check before building anything.
        let budget = 6 * n + 1024;
        let mut count = 0usize;
        for i in 0..n {
            cells(i, &mut |_, _| count += 1);
            if count > budget {
                return None;
            }
        }
        let mut pairs: Vec<(u32, u32)> = Vec::with_capacity(count);
        for i in 0..n {
            cells(i, &mut |cx, cy| {
                pairs.push(((cy * nx + cx) as u32, i as u32))
            });
        }
        let (start, items) = csr(nx * ny, &pairs);
        // Too dense somewhere: the adaptive partition will do better.
        if start.windows(2).any(|w| w[1] - w[0] > 256) {
            return None;
        }
        let mut leaves = Vec::with_capacity(nx * ny);
        for cy in 0..ny as i64 {
            for cx in 0..nx as i64 {
                let lo = Point::new(x0 + cx * s, y0 + cy * s);
                leaves.push((lo, Point::new(lo.x + s, lo.y + s)));
            }
        }
        Some(Grid {
            leaves,
            start,
            items,
        })
    }

    fn build_kd(segs: &[(Point, Point)], bboxes: &[Rect]) -> Grid {
        let mut g = Grid {
            leaves: Vec::new(),
            start: vec![0],
            items: Vec::new(),
        };
        let Some(bb) = bboxes.iter().copied().reduce(|a, b| a.union(&b)) else {
            g.leaves.push((Point::new(0, 0), Point::new(1, 1)));
            g.start.push(0);
            return g;
        };
        let lo = bb.min;
        let hi = Point::new(bb.max.x + 1, bb.max.y + 1);
        let all: Vec<u32> = (0..segs.len() as u32).collect();
        let mut stack: Vec<(Vec<u32>, Point, Point, u32)> = vec![(all, lo, hi, 0)];
        let mut sample: Vec<i64> = Vec::new();
        while let Some((items, lo, hi, depth)) = stack.pop() {
            if items.len() > LEAF_SIZE
                && depth < MAX_DEPTH
                && let Some((axis, at, left, right)) =
                    Self::split(segs, bboxes, &items, lo, hi, &mut sample)
            {
                let (lhi, rlo) = if axis == 0 {
                    (Point::new(at, hi.y), Point::new(at, lo.y))
                } else {
                    (Point::new(hi.x, at), Point::new(lo.x, at))
                };
                drop(items);
                stack.push((right, rlo, hi, depth + 1));
                stack.push((left, lo, lhi, depth + 1));
                continue;
            }
            g.items.extend_from_slice(&items);
            g.start.push(g.items.len() as u32);
            g.leaves.push((lo, hi));
        }
        g
    }

    /// Chooses a split of the region `[lo, hi)` and distributes the items, or `None` when no
    /// split makes progress.
    #[allow(clippy::type_complexity)]
    fn split(
        segs: &[(Point, Point)],
        bboxes: &[Rect],
        items: &[u32],
        lo: Point,
        hi: Point,
        sample: &mut Vec<i64>,
    ) -> Option<(u8, i64, Vec<u32>, Vec<u32>)> {
        let n = items.len();
        let stride = (n / 64).max(1);
        let mut best: Option<(usize, u8, i64)> = None;
        for axis in [0u8, 1u8] {
            let (l, h) = (coord(lo, axis), coord(hi, axis));
            if h - l < 2 {
                continue;
            }
            sample.clear();
            sample.extend(items.iter().step_by(stride).map(|&i| {
                let b = &bboxes[i as usize];
                let (a, c) = (coord(b.min, axis), coord(b.max, axis));
                a + (c - a) / 2
            }));
            let mid = sample.len() / 2;
            let mut at = *sample.select_nth_unstable(mid).1;
            if at <= l || at >= h {
                at = l + (h - l) / 2;
            }
            let at = at.clamp(l + 1, h - 1);
            // Straddlers in the sample (bounding box reaching both children).
            let straddle = items
                .iter()
                .step_by(stride)
                .filter(|&&i| {
                    let b = &bboxes[i as usize];
                    coord(b.min, axis) <= at && coord(b.max, axis) >= at - 1
                })
                .count();
            if best.is_none_or(|(s, _, _)| straddle < s) {
                best = Some((straddle, axis, at));
            }
        }
        let (_, axis, at) = best?;
        // Child rectangles grown by one unit (closed): left covers [lo - 1, at], right
        // [at - 1, hi] along the axis.
        let (lrect, rrect) = if axis == 0 {
            (
                Rect {
                    min: Point::new(lo.x - 1, lo.y - 1),
                    max: Point::new(at, hi.y),
                },
                Rect {
                    min: Point::new(at - 1, lo.y - 1),
                    max: hi,
                },
            )
        } else {
            (
                Rect {
                    min: Point::new(lo.x - 1, lo.y - 1),
                    max: Point::new(hi.x, at),
                },
                Rect {
                    min: Point::new(lo.x - 1, at - 1),
                    max: hi,
                },
            )
        };
        let mut left = Vec::with_capacity(n / 2 + 4);
        let mut right = Vec::with_capacity(n / 2 + 4);
        for &i in items {
            let b = &bboxes[i as usize];
            let (bmin, bmax) = (coord(b.min, axis), coord(b.max, axis));
            if bmax < at - 1 {
                left.push(i);
            } else if bmin > at {
                right.push(i);
            } else {
                let s = &segs[i as usize];
                if segment_meets_rect(s.0, s.1, &lrect) {
                    left.push(i);
                }
                if segment_meets_rect(s.0, s.1, &rrect) {
                    right.push(i);
                }
            }
        }
        // No progress: most items went to both sides.
        if (left.len() + right.len()) * 10 > n * 15 {
            return None;
        }
        Some((axis, at, left, right))
    }

    #[inline]
    fn n_cells(&self) -> usize {
        self.leaves.len()
    }

    /// `true` when `p` lies in leaf `c`'s half-open region.
    #[inline]
    fn in_leaf(&self, c: usize, p: Point) -> bool {
        let (lo, hi) = self.leaves[c];
        p.x >= lo.x && p.x < hi.x && p.y >= lo.y && p.y < hi.y
    }

    #[inline]
    fn items(&self, c: usize) -> &[u32] {
        &self.items[self.start[c] as usize..self.start[c + 1] as usize]
    }
}

/// Calls `f(i, j)` for every pair of segments in the cell whose bounding boxes overlap.
fn for_each_pair(
    segs: &[(Point, Point)],
    items: &[u32],
    bboxes: &[Rect],
    order: &mut Vec<(i128, i128, u32)>,
    mut f: impl FnMut(u32, u32),
) {
    pairs_dyn(segs, items, bboxes, order, &mut f)
}

fn pairs_dyn(
    segs: &[(Point, Point)],
    items: &[u32],
    bboxes: &[Rect],
    order: &mut Vec<(i128, i128, u32)>,
    f: &mut dyn FnMut(u32, u32),
) {
    if items.len() <= 24 {
        for (k, &i) in items.iter().enumerate() {
            let bi = &bboxes[i as usize];
            for &j in &items[k + 1..] {
                if bi.intersects(&bboxes[j as usize]) {
                    f(i, j);
                }
            }
        }
        return;
    }
    // A hub — an endpoint shared by many segments of the leaf — makes every pair of its
    // segments overlap in any projection. Two segments sharing an endpoint can only meet
    // elsewhere if they are collinear in the same direction, so hub pairs are found by
    // sorting the hub's segments by angle; the rest goes through the sweeps.
    if items.len() > 64
        && let Some((hub, h, rest)) = split_hub(segs, items)
    {
        let dir = |i: u32| {
            let (a, b) = segs[i as usize];
            if a == hub { sub(b, a) } else { sub(a, b) }
        };
        let mut hs = h;
        hs.sort_unstable_by(|&x, &y| crate::predicates::cmp_angle(dir(x), dir(y)).then(x.cmp(&y)));
        let mut k = 0;
        while k < hs.len() {
            let mut e = k + 1;
            while e < hs.len()
                && crate::predicates::cross(dir(hs[k]), dir(hs[e])) == 0
                && dot(dir(hs[k]), dir(hs[e])) > 0
            {
                e += 1;
            }
            for x in k..e {
                for y in x + 1..e {
                    f(hs[x], hs[y]);
                }
            }
            k = e;
        }
        pairs_dyn(segs, &rest, bboxes, order, f);
        bipartite_pairs(&hs, &rest, bboxes, f);
        return;
    }
    // Sweep-and-prune along the direction where the segments are thinnest: segments that
    // meet have overlapping projections on every direction.
    let d = Dir::best(segs, items);
    order.clear();
    order.extend(items.iter().map(|&i| {
        let (lo, hi) = d.range(&segs[i as usize]);
        (lo, hi, i)
    }));
    order.sort_unstable();
    for (k, &(_, hi, i)) in order.iter().enumerate() {
        let bi = &bboxes[i as usize];
        for &(lo2, _, j) in &order[k + 1..] {
            if lo2 > hi {
                break;
            }
            if bi.intersects(&bboxes[j as usize]) {
                f(i, j);
            }
        }
    }
}

/// The most frequent endpoint of `items` if it is shared by many of them: (hub, segments
/// incident to it, the others).
fn split_hub(segs: &[(Point, Point)], items: &[u32]) -> Option<(Point, Vec<u32>, Vec<u32>)> {
    let mut ends: Vec<Point> = items
        .iter()
        .flat_map(|&i| [segs[i as usize].0, segs[i as usize].1])
        .collect();
    ends.sort_unstable();
    let (mut best, mut best_n) = (ends[0], 0usize);
    let mut k = 0;
    while k < ends.len() {
        let mut e = k;
        while e < ends.len() && ends[e] == ends[k] {
            e += 1;
        }
        if e - k > best_n {
            (best, best_n) = (ends[k], e - k);
        }
        k = e;
    }
    if best_n < 16 || best_n * 4 < items.len() {
        return None;
    }
    let (h, rest): (Vec<u32>, Vec<u32>) = items
        .iter()
        .partition(|&&i| segs[i as usize].0 == best || segs[i as usize].1 == best);
    Some((best, h, rest))
}

/// Calls `f` for every pair across `a` and `b` whose bounding boxes overlap (sweep along x).
fn bipartite_pairs(a: &[u32], b: &[u32], bboxes: &[Rect], f: &mut dyn FnMut(u32, u32)) {
    let mut ev: Vec<(i64, bool, u32)> = a
        .iter()
        .map(|&i| (bboxes[i as usize].min.x, false, i))
        .collect();
    ev.extend(b.iter().map(|&i| (bboxes[i as usize].min.x, true, i)));
    ev.sort_unstable();
    let mut act: [Vec<u32>; 2] = [Vec::new(), Vec::new()];
    for (x, side, i) in ev {
        let other = &mut act[!side as usize];
        other.retain(|&j| bboxes[j as usize].max.x >= x);
        let bi = &bboxes[i as usize];
        for &j in other.iter() {
            if bi.intersects(&bboxes[j as usize]) {
                f(i, j);
            }
        }
        act[side as usize].push(i);
    }
}

/// Exact `floor(v + 1/2)` of `v = base + num * d / den` (`den > 0`).
#[inline]
fn round_exact(base: i64, num: i128, d: i128, den: i128) -> i64 {
    base + floor_div(2 * num * d + den, 2 * den) as i64
}

/// The proper crossing of `a-b` and `c-d`, rounded to the nearest integer point (ties up),
/// using floating point when it is provably correct and exact arithmetic otherwise.
#[inline]
fn rounded_crossing(a: Point, b: Point, c: Point, d: Point) -> Point {
    let o3 = orient(c, d, a);
    let o4 = orient(c, d, b);
    let (num, den) = if o3 - o4 < 0 {
        (-o3, o4 - o3)
    } else {
        (o3, o3 - o4)
    };
    let t = num as f64 / den as f64;
    let dx = (b.x - a.x) as i128;
    let dy = (b.y - a.y) as i128;
    // Absolute error of the float evaluation is below 2^-9 for |coords| <= 2^40.
    let round = |base: i64, dv: i128| -> i64 {
        let v = base as f64 + dv as f64 * t + 0.5;
        let f = libm::floor(v);
        let frac = v - f;
        if frac > 0.01 && frac < 0.99 {
            f as i64
        } else {
            round_exact(base, num, dv, den)
        }
    };
    Point::new(round(a.x, dx), round(a.y, dy))
}

/// Snap-rounds `segs` (each with distinct endpoints). Returns the fragments of every
/// segment, grouped by segment in input order and ordered from `a` to `b` within a segment.
pub(crate) fn snap_round(segs: &[(Point, Point)]) -> Vec<Frag> {
    let bboxes: Vec<Rect> = segs.iter().map(seg_bbox).collect();
    let grid = Grid::build(segs, &bboxes);
    let pp = pair_pass(segs, &bboxes, &grid);
    let mut st = Snap::new(segs, &bboxes, &grid, &pp);
    st.run();
    st.fragments(&pp.splits)
}

/// Output of the pair pass.
struct PairPass {
    /// Candidate pixels (endpoints and rounded crossings) per leaf (CSR), sorted by
    /// `(x, y)` without duplicates.
    pstart: Vec<u32>,
    pix: Vec<Point>,
    /// Crossing pixels (initially hot).
    crossing: Vec<bool>,
    /// Exact T-junctions: segment, squared distance from its start, point.
    splits: Csr2,
}

type Csr = (Vec<u32>, Vec<u32>);
type Csr2 = (Vec<u32>, Vec<(i128, Point)>);

/// Finds proper crossings (as rounded pixels) and exact T-junctions (a segment endpoint
/// on another segment's interior), and collects the candidate pixels of every leaf.
#[inline(never)]
fn pair_pass(segs: &[(Point, Point)], bboxes: &[Rect], grid: &Grid) -> PairPass {
    type Chunk = (Vec<(u32, (Point, bool))>, Vec<(u32, (i128, Point))>);
    let chunks: Vec<Chunk> = crate::par::map_ranges(grid.n_cells(), |range| {
        let mut order = Vec::new();
        let mut hot: Vec<(u32, (Point, bool))> = Vec::new();
        let mut splits: Vec<(u32, (i128, Point))> = Vec::new();
        for c in range {
            let items = grid.items(c);
            // Endpoints, each reported by the leaf containing it.
            for &i in items {
                let (a, b) = segs[i as usize];
                for p in [a, b] {
                    if grid.in_leaf(c, p) {
                        hot.push((c as u32, (p, false)));
                    }
                }
            }
            if items.len() < 2 {
                continue;
            }
            for_each_pair(segs, items, bboxes, &mut order, |i, j| {
                let (a, b) = segs[i as usize];
                let (p, q) = segs[j as usize];
                let o1 = orient(a, b, p).signum();
                let o2 = orient(a, b, q).signum();
                let o3 = orient(p, q, a).signum();
                let o4 = orient(p, q, b).signum();
                if o1 * o2 < 0 && o3 * o4 < 0 {
                    let x = rounded_crossing(a, b, p, q);
                    // Report each crossing once: in the leaf containing its pixel.
                    if grid.in_leaf(c, x) {
                        hot.push((c as u32, (x, true)));
                    }
                    return;
                }
                // T-junctions (each reported by the leaf containing the junction point).
                let mut tj = |s: u32, u: Point, v: Point, o: i32, x: Point| {
                    if o == 0
                        && x != u
                        && x != v
                        && grid.in_leaf(c, x)
                        && crate::predicates::on_segment(u, v, x)
                    {
                        splits.push((s, (crate::predicates::dist2(u, x), x)));
                    }
                };
                tj(i, a, b, o1 as i32, p);
                tj(i, a, b, o2 as i32, q);
                tj(j, p, q, o3 as i32, a);
                tj(j, p, q, o4 as i32, b);
            });
        }
        (hot, splits)
    });
    // Concatenate in leaf order (deterministic for any number of threads).
    let mut hot: Vec<(u32, (Point, bool))> =
        Vec::with_capacity(chunks.iter().map(|c| c.0.len()).sum());
    let mut splits: Vec<(u32, (i128, Point))> = Vec::new();
    for (h, s) in chunks {
        hot.extend(h);
        splits.extend(s);
    }
    let (pstart, mut pix) = csr(grid.n_cells(), &hot);
    drop(hot);
    // Sort and dedup within each leaf (merging flags), compacting in place.
    let mut out_start = vec![0u32; grid.n_cells() + 1];
    let mut w = 0usize;
    for c in 0..grid.n_cells() {
        let (s0, s1) = (pstart[c] as usize, pstart[c + 1] as usize);
        pix[s0..s1].sort_unstable();
        for k in s0..s1 {
            let (p, f) = pix[k];
            if w > out_start[c] as usize && pix[w - 1].0 == p {
                pix[w - 1].1 |= f;
            } else {
                pix[w] = (p, f);
                w += 1;
            }
        }
        out_start[c + 1] = w as u32;
    }
    pix.truncate(w);
    PairPass {
        pstart: out_start,
        crossing: pix.iter().map(|x| x.1).collect(),
        pix: pix.into_iter().map(|x| x.0).collect(),
        splits: csr(segs.len(), &splits),
    }
}

/// How a segment relates to a candidate pixel.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Rel {
    /// Too far to matter.
    Far,
    /// The segment meets the pixel.
    Meets,
    /// The pixel centre is within `sqrt(2)/2` of the segment's line but the segment
    /// misses the pixel: the centre could still lie on a rerouted fragment.
    Near,
}

#[inline]
fn relation(a: Point, b: Point, bb: &Rect, p: Point) -> Rel {
    if p.x < bb.min.x - 1 || p.x > bb.max.x + 1 || p.y < bb.min.y - 1 || p.y > bb.max.y + 1 {
        return Rel::Far;
    }
    // Meeting the pixel, or lying on a fragment, requires the centre to be within
    // sqrt(2)/2 of the segment's line: |cross| <= 0.71 * len. The float cross product is
    // accurate to far better than the slack used here.
    let dx = (b.x - a.x) as f64;
    let dy = (b.y - a.y) as f64;
    let len2 = dx * dx + dy * dy;
    let px = (p.x - a.x) as f64;
    let py = (p.y - a.y) as f64;
    let o = dx * py - dy * px;
    if o * o > 0.55 * len2 {
        return Rel::Far;
    }
    // Fast accept: within 1/2 of the line and well inside the segment's span means the
    // segment crosses the pixel's inscribed disk.
    let t = dx * px + dy * py;
    let l = len2.sqrt();
    if (o * o < 0.24 * len2 && t > l * 0.75 && t < len2 - l * 0.75)
        || segment_pixel_entry(a, b, p).is_some()
    {
        Rel::Meets
    } else {
        Rel::Near
    }
}

/// Snap-rounding state: which candidate pixels are hot and which segments are rerouted.
///
/// Crossing pixels are hot. A segment meeting a hot pixel (other than at its own endpoints)
/// is rerouted ("affected"), and every candidate pixel an affected segment meets, including
/// its own endpoints, becomes hot in turn — until nothing changes. Segments that are not
/// affected keep their exact original position, so input without crossings is returned
/// unchanged (only split at vertices lying exactly on segments). Within the affected part,
/// this is Hobby's snap rounding with all relevant hot pixels, which preserves topology.
///
/// Work is done leaf by leaf, lazily: a leaf's segment/pixel relations are computed the
/// first time one of its pixels becomes hot or one of its segments becomes affected, so the
/// cost follows the regions where rounding actually happens.
struct Snap<'a> {
    segs: &'a [(Point, Point)],
    bboxes: &'a [Rect],
    grid: &'a Grid,
    pstart: &'a [u32],
    pix: &'a [Point],
    /// Leaf of every candidate pixel.
    pleaf: Vec<u32>,
    /// Leaves of every segment (CSR).
    sleaves: Csr,
    processed: Vec<bool>,
    active: Vec<bool>,
    affected: Vec<bool>,
    /// Relations found so far, as linked lists: `rel[i] = (segment, pixel, next for the
    /// segment, next for the pixel)`; `seg_head` / `pix_head` start the lists.
    rel: Vec<(u32, u32, u32, u32)>,
    seg_head: Vec<u32>,
    pix_head: Vec<u32>,
    /// Pixels met by segments (including own endpoints, marked by `own`), for fragments.
    near: Vec<(u32, u32)>,
    leaf_queue: Vec<u32>,
    px_queue: Vec<u32>,
    seg_queue: Vec<u32>,
    by_d: Vec<(i128, u32)>,
}

const NIL: u32 = u32::MAX;

impl<'a> Snap<'a> {
    fn new(
        segs: &'a [(Point, Point)],
        bboxes: &'a [Rect],
        grid: &'a Grid,
        pp: &'a PairPass,
    ) -> Self {
        let mut pleaf = vec![0u32; pp.pix.len()];
        for c in 0..grid.n_cells() {
            for k in pp.pstart[c]..pp.pstart[c + 1] {
                pleaf[k as usize] = c as u32;
            }
        }
        let mut sl: Vec<(u32, u32)> = Vec::with_capacity(grid.items.len());
        for c in 0..grid.n_cells() {
            for &s in grid.items(c) {
                sl.push((s, c as u32));
            }
        }
        let mut leaf_queue: Vec<u32> = (0..pp.pix.len())
            .filter(|&k| pp.crossing[k])
            .map(|k| pleaf[k])
            .collect();
        leaf_queue.dedup();
        Snap {
            segs,
            bboxes,
            grid,
            pstart: &pp.pstart,
            pix: &pp.pix,
            pleaf,
            sleaves: csr(segs.len(), &sl),
            processed: vec![false; grid.n_cells()],
            active: pp.crossing.clone(),
            affected: vec![false; segs.len()],
            rel: Vec::new(),
            seg_head: vec![NIL; segs.len()],
            pix_head: vec![NIL; pp.pix.len()],
            near: Vec::new(),
            leaf_queue,
            px_queue: Vec::new(),
            seg_queue: Vec::new(),
            by_d: Vec::new(),
        }
    }

    /// Segment/pixel relations of leaf `c` (pure; `by_d` is scratch space).
    fn leaf_relations(
        &self,
        c: usize,
        by_d: &mut Vec<(i128, u32)>,
        found: &mut Vec<(u32, u32, Rel)>,
    ) {
        let base = self.pstart[c] as usize;
        let cp = &self.pix[base..self.pstart[c + 1] as usize];
        if cp.is_empty() {
            return;
        }
        let items = self.grid.items(c);
        // For leaves with many pixels, sort them along the direction where the segments
        // are thinnest (pixels are already sorted by x).
        let dir = if cp.len() > 16 {
            Dir::best(self.segs, items)
        } else {
            Dir::X
        };
        by_d.clear();
        if dir != Dir::X {
            by_d.extend(cp.iter().enumerate().map(|(k, &p)| (dir.proj(p), k as u32)));
            by_d.sort_unstable();
        }
        for &s in items {
            let (a, b) = self.segs[s as usize];
            let bb = &self.bboxes[s as usize];
            let mut test = |k: usize| {
                let p = cp[k];
                let id = (base + k) as u32;
                if p == a || p == b {
                    // Own endpoint: relevant only once the segment is affected.
                    found.push((s, id, Rel::Far));
                    return;
                }
                let r = relation(a, b, bb, p);
                if r != Rel::Far {
                    found.push((s, id, r));
                }
            };
            if cp.len() <= 16 {
                (0..cp.len()).for_each(&mut test);
            } else if dir == Dir::X {
                let from = cp.partition_point(|p| p.x < bb.min.x - 1);
                let to = from + cp[from..].partition_point(|p| p.x <= bb.max.x + 1);
                (from..to).for_each(&mut test);
            } else {
                // Pixels within distance 1 of the segment project within the margin of its
                // projected range.
                let (lo, hi) = dir.range(&(a, b));
                let m = dir.margin();
                let from = by_d.partition_point(|x| x.0 < lo - m);
                let to = from + by_d[from..].partition_point(|x| x.0 <= hi + m);
                by_d[from..to].iter().for_each(|&(_, k)| test(k as usize));
            }
        }
    }

    /// Computes all segment/pixel relations of leaf `c` and applies them.
    fn process_leaf(&mut self, c: usize) {
        if self.processed[c] {
            return;
        }
        let mut by_d = core::mem::take(&mut self.by_d);
        let mut found = Vec::new();
        self.leaf_relations(c, &mut by_d, &mut found);
        self.by_d = by_d;
        self.apply(c, found);
    }

    /// Records the relations of leaf `c`, propagating hot pixels and affected segments.
    fn apply(&mut self, c: usize, found: Vec<(u32, u32, Rel)>) {
        self.processed[c] = true;
        for (s, k, r) in found {
            match r {
                Rel::Near => self.near.push((s, k)),
                // Meets, or own endpoint (`Far` marker): linked into both lists; own
                // endpoints are distinguished when used.
                _ => {
                    let i = self.rel.len() as u32;
                    self.rel
                        .push((s, k, self.seg_head[s as usize], self.pix_head[k as usize]));
                    self.seg_head[s as usize] = i;
                    self.pix_head[k as usize] = i;
                    let own = {
                        let (a, b) = self.segs[s as usize];
                        self.pix[k as usize] == a || self.pix[k as usize] == b
                    };
                    if self.affected[s as usize] {
                        self.activate(k);
                    } else if !own && self.active[k as usize] {
                        self.affect(s);
                    }
                }
            }
        }
    }

    /// When crossings are dense almost every leaf ends up processed: compute all leaf
    /// relations up front (in parallel with the `rayon` feature) and apply them in leaf
    /// order. The fixpoint reached is the same as with lazy processing.
    fn precompute_all(&mut self) {
        let n = self.grid.n_cells();
        let this = &*self;
        // (first leaf, relations of the chunk's leaves, end offset of each leaf's relations)
        type Chunk = (usize, Vec<(u32, u32, Rel)>, Vec<u32>);
        let chunks: Vec<Chunk> = crate::par::map_ranges(n, |range| {
            let mut by_d = Vec::new();
            let mut found = Vec::new();
            let mut ends = Vec::with_capacity(range.len());
            let start = range.start;
            for c in range {
                this.leaf_relations(c, &mut by_d, &mut found);
                ends.push(found.len() as u32);
            }
            (start, found, ends)
        });
        for (start, found, ends) in chunks {
            let mut it = found.into_iter();
            let mut prev = 0u32;
            for (i, &e) in ends.iter().enumerate() {
                let part: Vec<(u32, u32, Rel)> = it.by_ref().take((e - prev) as usize).collect();
                prev = e;
                self.apply(start + i, part);
            }
        }
    }

    fn activate(&mut self, k: u32) {
        if !self.active[k as usize] {
            self.active[k as usize] = true;
            self.px_queue.push(k);
        }
    }

    fn affect(&mut self, s: u32) {
        if !self.affected[s as usize] {
            self.affected[s as usize] = true;
            self.seg_queue.push(s);
        }
    }

    fn run(&mut self) {
        // Crossing pixels are hot from the start: process their leaves, then propagate.
        let mut hot = 0usize;
        for k in 0..self.active.len() {
            if self.active[k] {
                self.px_queue.push(k as u32);
                hot += 1;
            }
        }
        if hot * 16 > self.segs.len() {
            self.precompute_all();
        }
        loop {
            if let Some(c) = self.leaf_queue.pop() {
                self.process_leaf(c as usize);
            } else if let Some(k) = self.px_queue.pop() {
                let c = self.pleaf[k as usize] as usize;
                if !self.processed[c] {
                    // Processing the leaf will apply this pixel's relations.
                    self.process_leaf(c);
                }
                let p = self.pix[k as usize];
                let mut i = self.pix_head[k as usize];
                while i != NIL {
                    let (s, _, _, next) = self.rel[i as usize];
                    let (a, b) = self.segs[s as usize];
                    if p != a && p != b {
                        self.affect(s);
                    }
                    i = next;
                }
            } else if let Some(s) = self.seg_queue.pop() {
                // All leaves of an affected segment must be known, then everything it meets
                // (and its own endpoints) becomes hot.
                let (ls, le) = (
                    self.sleaves.0[s as usize] as usize,
                    self.sleaves.0[s as usize + 1] as usize,
                );
                for li in ls..le {
                    let c = self.sleaves.1[li] as usize;
                    self.process_leaf(c);
                }
                let mut i = self.seg_head[s as usize];
                while i != NIL {
                    let (_, k, next, _) = self.rel[i as usize];
                    self.activate(k);
                    i = next;
                }
            } else {
                break;
            }
        }
    }

    /// Builds the fragments of every segment.
    fn fragments(mut self, splits: &Csr2) -> Vec<Frag> {
        let n = self.segs.len();
        let mut hit_pairs: Vec<(u32, u32)> = Vec::new();
        for s in 0..n {
            if !self.affected[s] {
                continue;
            }
            let (a, b) = self.segs[s];
            let mut i = self.seg_head[s];
            while i != NIL {
                let (_, k, next, _) = self.rel[i as usize];
                let p = self.pix[k as usize];
                if p != a && p != b {
                    hit_pairs.push((s as u32, k));
                }
                i = next;
            }
        }
        let (hstart, hits) = csr(n, &hit_pairs);
        let (nstart, near) = csr(n, &core::mem::take(&mut self.near));
        let (sstart, sp) = (&splits.0, &splits.1);
        let mut out: Vec<Frag> = Vec::with_capacity(n + hits.len() + sp.len());
        let mut poly: Vec<Point> = Vec::new();
        let mut ord: Vec<(i128, Point)> = Vec::new();
        let mut ins: Vec<(usize, i128, Point)> = Vec::new();
        for (si, &(a, b)) in self.segs.iter().enumerate() {
            let src = si as u32;
            if !self.affected[si] {
                // Unaffected: exact position, split at T-junctions only.
                let s = &sp[sstart[si] as usize..sstart[si + 1] as usize];
                if s.is_empty() {
                    out.push(Frag { a, b, src });
                    continue;
                }
                ord.clear();
                ord.extend_from_slice(s);
                ord.sort_unstable();
                let mut cur = a;
                for &(_, x) in &ord {
                    if x != cur {
                        out.push(Frag { a: cur, b: x, src });
                        cur = x;
                    }
                }
                out.push(Frag { a: cur, b, src });
                continue;
            }
            // Hot pixels are met in the order of their centres' projections on the segment
            // direction (pixel rows and columns are traversed monotonically).
            let h = &hits[hstart[si] as usize..hstart[si + 1] as usize];
            let d = sub(b, a);
            ord.clear();
            ord.extend(
                h.iter()
                    .map(|&k| self.pix[k as usize])
                    .map(|p| (dot(sub(p, a), d), p)),
            );
            ord.sort_unstable();
            ord.dedup();
            poly.clear();
            poly.push(a);
            poly.extend(ord.iter().map(|x| x.1));
            poly.push(b);
            // T-junctions: candidate centres lying on a rerouted fragment's interior.
            ins.clear();
            for &k in &near[nstart[si] as usize..nstart[si + 1] as usize] {
                let c = self.pix[k as usize];
                for k in 0..poly.len() - 1 {
                    if in_segment_interior(poly[k], poly[k + 1], c) {
                        ins.push((k, crate::predicates::dist2(poly[k], c), c));
                        break;
                    }
                }
            }
            ins.sort_unstable();
            ins.dedup();
            let mut j = 0;
            for k in 0..poly.len() - 1 {
                let mut cur = poly[k];
                while j < ins.len() && ins[j].0 == k {
                    if ins[j].2 != cur {
                        out.push(Frag {
                            a: cur,
                            b: ins[j].2,
                            src,
                        });
                        cur = ins[j].2;
                    }
                    j += 1;
                }
                if poly[k + 1] != cur {
                    out.push(Frag {
                        a: cur,
                        b: poly[k + 1],
                        src,
                    });
                }
            }
        }
        out
    }
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
    let n = segs.len();
    let bboxes: Vec<Rect> = segs.iter().map(seg_bbox).collect();
    let grid = Grid::build(segs, &bboxes);
    let mut order = Vec::new();
    let mut crossing: Option<Crossing> = None;
    let mut splits: Vec<(u32, (i128, Point))> = Vec::new();
    for c in 0..grid.n_cells() {
        let items = grid.items(c);
        if items.len() < 2 {
            continue;
        }
        for_each_pair(segs, items, &bboxes, &mut order, |i, j| {
            let (a, b) = segs[i as usize];
            let (p, q) = segs[j as usize];
            if segments_cross_properly(a, b, p, q) {
                let cr = Crossing {
                    i: i.min(j),
                    j: i.max(j),
                };
                if crossing.is_none_or(|c0| cr < c0) {
                    crossing = Some(cr);
                }
                return;
            }
            for (s, (u, v), pts) in [(i, (a, b), [p, q]), (j, (p, q), [a, b])] {
                for x in pts {
                    if grid.in_leaf(c, x) && in_segment_interior(u, v, x) {
                        splits.push((s, (crate::predicates::dist2(u, x), x)));
                    }
                }
            }
        });
    }
    if let Some(c) = crossing {
        return Err(c);
    }
    let (sstart, mut sp) = csr(n, &splits);
    drop(splits);
    let mut out = Vec::with_capacity(n + sp.len());
    for (si, &(a, b)) in segs.iter().enumerate() {
        let s = &mut sp[sstart[si] as usize..sstart[si + 1] as usize];
        s.sort_unstable();
        let si = si as u32;
        let mut cur = a;
        for &(_, x) in s.iter() {
            if x != cur {
                out.push(Frag {
                    a: cur,
                    b: x,
                    src: si,
                });
                cur = x;
            }
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

    #[test]
    fn random_hub_noded() {
        // Many segments sharing one endpoint, some collinear (overlapping), plus others.
        let mut s: u64 = 77;
        let mut rnd = |m: i64| {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((s >> 33) as i64).rem_euclid(m)
        };
        for _ in 0..60 {
            let hub = p(rnd(50), rnd(50));
            let mut segs: Vec<(Point, Point)> = (0..120)
                .map(|_| {
                    let (dx, dy) = (rnd(9) - 4, rnd(9) - 4);
                    let k = 1 + rnd(12);
                    (hub, p(hub.x + dx * k, hub.y + dy * k))
                })
                .filter(|(a, b)| a != b)
                .collect();
            for _ in 0..30 {
                segs.push((
                    p(rnd(100) - 25, rnd(100) - 25),
                    p(rnd(100) - 25, rnd(100) - 25),
                ));
            }
            segs.retain(|(a, b)| a != b);
            assert_noded(&snap_round(&segs));
            if let Ok(f) = node_exact(&segs) {
                assert_noded(&f);
            }
        }
    }

    #[test]
    fn random_axis_parallel_noded() {
        let mut s: u64 = 0xdead_beef;
        let mut rnd = |m: i64| {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((s >> 33) as i64).rem_euclid(m)
        };
        for it in 0..400 {
            let n = 2 + rnd(80) as usize;
            let range = 2 + rnd(if it % 2 == 0 { 20 } else { 2000 });
            let segs: Vec<(Point, Point)> = (0..n)
                .map(|_| {
                    let a = p(rnd(range), rnd(range));
                    match rnd(3) {
                        0 => (a, p(a.x, rnd(range))),
                        1 => (a, p(rnd(range), a.y)),
                        _ => (a, p(rnd(range), rnd(range))),
                    }
                })
                .filter(|(a, b)| a != b)
                .collect();
            for f in [snap_round(&segs), node_exact(&segs).unwrap_or_default()] {
                assert_noded(&f);
            }
        }
    }
}

#[cfg(test)]
mod regression {
    use super::*;
    #[test]
    fn degenerate_parallelograms_noded() {
        // Long vertical segments were once registered in the wrong grid cells.
        let p = Point::new;
        let sq =
            |x0: i64, y0: i64, x1: i64, y1: i64| vec![p(x0, y0), p(x1, y0), p(x1, y1), p(x0, y1)];
        let a_rings = vec![sq(0, 0, 10, 10), vec![p(4, 4), p(4, 6), p(6, 6), p(6, 4)]];
        let b = sq(0, 0, 1, 1);
        let mut rings: Vec<Vec<Point>> = Vec::new();
        for r in &a_rings {
            rings.push(r.clone());
        }
        for ra in &a_rings {
            let n = ra.len();
            for i in 0..n {
                let (p0, p1) = (ra[i], ra[(i + 1) % n]);
                rings.push(b.iter().map(|q| p(q.x + p0.x, q.y + p0.y)).collect());
                for j in 0..4 {
                    let (q0, q1) = (b[j], b[(j + 1) % 4]);
                    rings.push(vec![
                        p(p0.x + q0.x, p0.y + q0.y),
                        p(p1.x + q0.x, p1.y + q0.y),
                        p(p1.x + q1.x, p1.y + q1.y),
                        p(p0.x + q1.x, p0.y + q1.y),
                    ]);
                }
            }
        }
        let mut segs = Vec::new();
        for r in &rings {
            for i in 0..r.len() {
                let (a, b) = (r[i], r[(i + 1) % r.len()]);
                if a != b {
                    segs.push((a, b));
                }
            }
        }
        let frags = snap_round(&segs);
        for (i, f) in frags.iter().enumerate() {
            for g in &frags[i + 1..] {
                for v in [g.a, g.b] {
                    if in_segment_interior(f.a, f.b, v) {
                        panic!(
                            "{v:?} on {f:?} (seg {:?}) from {g:?} (seg {:?})",
                            segs[f.src as usize], segs[g.src as usize]
                        );
                    }
                }
                assert!(
                    !segments_cross_properly(f.a, f.b, g.a, g.b),
                    "{f:?} x {g:?}"
                );
            }
        }
    }
}
