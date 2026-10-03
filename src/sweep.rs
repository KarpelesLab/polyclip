//! Plane sweep over a fully noded arrangement.
//!
//! The sweep line moves in lexicographic `(x, y)` order (a vertical line sheared
//! infinitesimally so that vertical edges behave like edges leaning right). Because the
//! edges never cross and meet only at shared endpoints, their relative order in the status
//! never changes while they are active; at every vertex the edges ending there form a
//! contiguous run of the status, which is replaced by the edges starting there.

use crate::geom::Point;
use crate::predicates::{cmp_dir_halfplane, orient, sub};
use core::cmp::Ordering;

/// Orders edges `(lo, hi)` (with `lo < hi`) for [`sweep`]: by `lo`, then bottom to top.
#[inline]
pub(crate) fn cmp_sweep_edges(a: (Point, Point), b: (Point, Point)) -> Ordering {
    a.0.cmp(&b.0)
        .then_with(|| cmp_dir_halfplane(sub(a.1, a.0), sub(b.1, b.0)))
}

/// Sweep status: the active edges bottom to top, split into blocks of bounded size so that
/// insertions and removals stay cheap however wide the status grows. `last[b]` mirrors the
/// last edge of block `b` in a dense array for the top-level binary search.
struct Status {
    blocks: Vec<Vec<u32>>,
    last: Vec<u32>,
}

const BLOCK: usize = 256;

impl Status {
    fn new() -> Self {
        Status {
            blocks: Vec::new(),
            last: Vec::new(),
        }
    }

    /// Position `(block, index)` of the first edge for which `below` is false.
    #[inline]
    fn find(&self, below: impl Fn(u32) -> bool) -> (usize, usize) {
        let b = self.last.partition_point(|&e| below(e));
        if b == self.blocks.len() {
            return match self.blocks.last() {
                Some(l) => (b - 1, l.len()),
                None => (0, 0),
            };
        }
        (b, self.blocks[b].partition_point(|&e| below(e)))
    }

    /// The edge just before position `(b, i)`.
    #[inline]
    fn before(&self, b: usize, i: usize) -> Option<u32> {
        if i > 0 {
            Some(self.blocks[b][i - 1])
        } else if b > 0 {
            Some(self.last[b - 1])
        } else {
            None
        }
    }

    /// Removes `remove` edges at `(b, i)` and inserts `insert` there.
    fn splice(
        &mut self,
        b: usize,
        i: usize,
        mut remove: usize,
        insert: impl ExactSizeIterator<Item = u32>,
    ) {
        if self.blocks.is_empty() {
            if insert.len() == 0 {
                return;
            }
            self.blocks.push(Vec::with_capacity(2 * BLOCK));
            self.last.push(0);
        }
        // Remove (possibly spilling into following blocks).
        let mut bb = b;
        let mut ii = i;
        while remove > 0 && bb < self.blocks.len() {
            let blk = &mut self.blocks[bb];
            let k = remove.min(blk.len() - ii);
            blk.drain(ii..ii + k);
            remove -= k;
            bb += 1;
            ii = 0;
        }
        let blk = &mut self.blocks[b];
        let i = i.min(blk.len());
        blk.splice(i..i, insert);
        if blk.len() > 2 * BLOCK {
            let tail = blk.split_off(BLOCK);
            self.blocks.insert(b + 1, tail);
            self.last.insert(b + 1, 0);
        }
        // Refresh `last` for touched blocks and drop empty ones.
        let end = bb.max(b + 2).min(self.blocks.len());
        let mut k = b;
        let mut end = end;
        while k < end {
            if self.blocks[k].is_empty() {
                self.blocks.remove(k);
                self.last.remove(k);
                end -= 1;
            } else {
                self.last[k] = *self.blocks[k].last().unwrap();
                k += 1;
            }
        }
    }
}

/// Runs the sweep. `edges[i] = (lo, hi)` with `lo < hi`, sorted by [`cmp_sweep_edges`],
/// pairwise non-crossing, with no vertex on another edge's interior.
///
/// `on_insert(e, below)` is called for every edge, in sweep order and bottom to top at each
/// vertex, with the edge immediately below it at the moment it is inserted.
pub(crate) fn sweep(edges: &[(Point, Point)], mut on_insert: impl FnMut(u32, Option<u32>)) {
    let n = edges.len();
    let mut his: Vec<Point> = edges.iter().map(|e| e.1).collect();
    his.sort_unstable();
    let mut status = Status::new();
    let mut si = 0usize;
    let mut hi = 0usize;
    loop {
        let v = match (edges.get(si), his.get(hi)) {
            (Some(e), Some(&h)) => e.0.min(h),
            (Some(e), None) => e.0,
            (None, Some(&h)) => h,
            (None, None) => break,
        };
        let mut ending = 0usize;
        while hi < n && his[hi] == v {
            ending += 1;
            hi += 1;
        }
        let s0 = si;
        while si < n && edges[si].0 == v {
            si += 1;
        }
        let (b, i) = status.find(|e| {
            let (lo, h) = edges[e as usize];
            h != v && orient(lo, h, v) > 0
        });
        let below = status.before(b, i);
        status.splice(b, i, ending, (s0..si).map(|e| e as u32));
        for k in s0..si {
            let bl = if k == s0 { below } else { Some((k - 1) as u32) };
            on_insert(k as u32, bl);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nested_squares() {
        let p = Point::new;
        // Outer square 0..10 and inner square 2..8, edges as (lo, hi).
        let mut e = vec![
            (p(0, 0), p(10, 0)),
            (p(0, 0), p(0, 10)),
            (p(0, 10), p(10, 10)),
            (p(10, 0), p(10, 10)),
            (p(2, 2), p(8, 2)),
            (p(2, 2), p(2, 8)),
            (p(2, 8), p(8, 8)),
            (p(8, 2), p(8, 8)),
        ];
        e.sort_by(|a, b| cmp_sweep_edges(*a, *b));
        let mut below = vec![None; e.len()];
        sweep(&e, |i, b| below[i as usize] = b.map(|b| e[b as usize]));
        let idx = |x: (Point, Point)| e.iter().position(|&y| y == x).unwrap();
        assert_eq!(below[idx((p(0, 0), p(10, 0)))], None);
        assert_eq!(below[idx((p(0, 0), p(0, 10)))], Some((p(0, 0), p(10, 0))));
        assert_eq!(below[idx((p(2, 2), p(8, 2)))], Some((p(0, 0), p(10, 0))));
        assert_eq!(below[idx((p(2, 8), p(8, 8)))], Some((p(2, 2), p(8, 2))));
        assert_eq!(below[idx((p(10, 0), p(10, 10)))], None);
    }
}
