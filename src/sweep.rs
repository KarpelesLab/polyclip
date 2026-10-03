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

/// Runs the sweep. `edges[i] = (lo, hi)` with `lo < hi`, sorted by [`cmp_sweep_edges`],
/// pairwise non-crossing, with no vertex on another edge's interior.
///
/// `on_insert(e, below)` is called for every edge, in sweep order and bottom to top at each
/// vertex, with the edge immediately below it at the moment it is inserted.
pub(crate) fn sweep(edges: &[(Point, Point)], mut on_insert: impl FnMut(u32, Option<u32>)) {
    let n = edges.len();
    let mut his: Vec<Point> = edges.iter().map(|e| e.1).collect();
    his.sort_unstable();
    let mut status: Vec<u32> = Vec::new();
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
        let pos = status.partition_point(|&e| {
            let (lo, h) = edges[e as usize];
            h != v && orient(lo, h, v) > 0
        });
        let end = (pos + ending).min(status.len());
        debug_assert!(
            status[pos..end].iter().all(|&e| edges[e as usize].1 == v),
            "sweep: non-noded input"
        );
        status.splice(pos..end, (s0..si).map(|e| e as u32));
        for k in s0..si {
            let below = if k == s0 {
                if pos > 0 { Some(status[pos - 1]) } else { None }
            } else {
                Some((k - 1) as u32)
            };
            on_insert(k as u32, below);
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
