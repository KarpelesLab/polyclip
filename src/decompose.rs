//! Trapezoidal (vertical) decomposition.

use crate::boolean::{FillRule, RingSource, union_all};
use crate::error::Result;
use crate::geom::{Point, PointF};
use core::cmp::Ordering;

/// A trapezoid of a vertical decomposition: the region with `x0 <= x <= x1` between the
/// lines supporting the `bottom` and `top` edges. Both edges span `[x0, x1]` and are given
/// left to right, so the representation is exact (corners are generally not integer
/// points).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Trapezoid {
    /// Left boundary.
    pub x0: i64,
    /// Right boundary (`> x0`).
    pub x1: i64,
    /// Lower boundary edge (left endpoint first).
    pub bottom: (Point, Point),
    /// Upper boundary edge (left endpoint first).
    pub top: (Point, Point),
}

fn y_at(e: (Point, Point), x: i64) -> f64 {
    let (a, b) = e;
    a.y as f64 + (x - a.x) as f64 * (b.y - a.y) as f64 / (b.x - a.x) as f64
}

impl Trapezoid {
    /// The four corners (lower-left, lower-right, upper-right, upper-left) as `f64`.
    pub fn corners(&self) -> [PointF; 4] {
        [
            PointF::new(self.x0 as f64, y_at(self.bottom, self.x0)),
            PointF::new(self.x1 as f64, y_at(self.bottom, self.x1)),
            PointF::new(self.x1 as f64, y_at(self.top, self.x1)),
            PointF::new(self.x0 as f64, y_at(self.top, self.x0)),
        ]
    }

    /// Area (as `f64`).
    pub fn area(&self) -> f64 {
        let w = (self.x1 - self.x0) as f64;
        let h0 = y_at(self.top, self.x0) - y_at(self.bottom, self.x0);
        let h1 = y_at(self.top, self.x1) - y_at(self.bottom, self.x1);
        w * (h0 + h1) / 2.0
    }
}

/// Compares edges `a` and `b` (left to right, non-vertical, non-crossing, both spanning the
/// slab `[xl, xr]`) by height at the slab's middle, exactly.
fn cmp_in_slab(a: (Point, Point), b: (Point, Point), xl: i64, xr: i64) -> Ordering {
    // y(x) = a0.y + (x - a0.x) * dy / dx at x = (xl + xr) / 2; compare doubled values:
    // (2 a0.y dx + (xl + xr - 2 a0.x) dy) / dx.
    let num = |e: (Point, Point)| -> (i128, i128) {
        let dx = (e.1.x - e.0.x) as i128;
        let dy = (e.1.y - e.0.y) as i128;
        (
            2 * e.0.y as i128 * dx + (xl as i128 + xr as i128 - 2 * e.0.x as i128) * dy,
            dx,
        )
    };
    let (na, da) = num(a);
    let (nb, db) = num(b);
    (na * db).cmp(&(nb * da))
}

/// Vertical decomposition of the region of `input` (non-zero fill rule) into trapezoids.
///
/// The region is first normalized with a union. Vertical walls are extended from every
/// vertex up and down to the nearest boundary edge, inside the region only; maximal pieces
/// between two walls with the same bottom and top edges form one trapezoid. The result is
/// sorted by `(x0, bottom, top)` and covers the region exactly with disjoint interiors.
/// Useful for free-space representation in routers (trapezoids are convex and adjacent
/// along vertical walls).
///
/// ```
/// use polyclip::{trapezoids, Polygon, Ring};
/// let l = Ring::from([(0, 0), (20, 0), (20, 10), (10, 10), (10, 20), (0, 20)]);
/// let t = trapezoids(&l).unwrap();
/// assert_eq!(t.len(), 2);
/// let area: f64 = t.iter().map(|t| t.area()).sum();
/// assert_eq!(area, 300.0);
/// ```
pub fn trapezoids(input: &(impl RingSource + ?Sized)) -> Result<Vec<Trapezoid>> {
    let polys = union_all(input, FillRule::NonZero)?;
    // Non-vertical edges, left to right.
    let mut edges: Vec<(Point, Point)> = Vec::new();
    for p in &polys {
        for r in p.rings() {
            for (a, b) in r.edges() {
                if a.x != b.x {
                    edges.push(if a.x < b.x { (a, b) } else { (b, a) });
                }
            }
        }
    }
    let mut xs: Vec<i64> = edges.iter().flat_map(|e| [e.0.x, e.1.x]).collect();
    xs.sort_unstable();
    xs.dedup();
    edges.sort_unstable_by_key(|e| (e.0.x, e.0.y, e.1.x, e.1.y));
    let mut out: Vec<Trapezoid> = Vec::new();
    // Open trapezoids: (bottom edge index, top edge index, start x).
    let mut open: Vec<(u32, u32, i64)> = Vec::new();
    let mut status: Vec<u32> = Vec::new();
    let mut next_edge = 0usize;
    for w in xs.windows(2) {
        let (xl, xr) = (w[0], w[1]);
        // Active edges for slab [xl, xr]: drop those ending at or before xl, add new ones.
        status.retain(|&e| edges[e as usize].1.x > xl);
        while next_edge < edges.len() && edges[next_edge].0.x <= xl {
            if edges[next_edge].1.x > xl {
                status.push(next_edge as u32);
            }
            next_edge += 1;
        }
        status.sort_by(|&a, &b| {
            cmp_in_slab(edges[a as usize], edges[b as usize], xl, xr).then(a.cmp(&b))
        });
        // Interior gaps are between status[2k] and status[2k + 1].
        let mut new_open: Vec<(u32, u32, i64)> = Vec::with_capacity(status.len() / 2);
        let mut k = 0;
        while k + 1 < status.len() {
            let (b, t) = (status[k], status[k + 1]);
            // Continue an open trapezoid with the same bounding edges, else start one.
            let start = match open.iter().position(|&(ob, ot, _)| ob == b && ot == t) {
                Some(i) => open.swap_remove(i).2,
                None => xl,
            };
            new_open.push((b, t, start));
            k += 2;
        }
        for (b, t, start) in open.drain(..) {
            out.push(Trapezoid {
                x0: start,
                x1: xl,
                bottom: edges[b as usize],
                top: edges[t as usize],
            });
        }
        open = new_open;
    }
    if let Some(&xlast) = xs.last() {
        for (b, t, start) in open.drain(..) {
            out.push(Trapezoid {
                x0: start,
                x1: xlast,
                bottom: edges[b as usize],
                top: edges[t as usize],
            });
        }
    }
    out.sort_unstable_by_key(|t| (t.x0, t.bottom, t.top, t.x1));
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::{Polygon, Ring};

    fn sq(x0: i64, y0: i64, x1: i64, y1: i64) -> Ring {
        Ring::from([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])
    }

    #[test]
    fn square_with_hole() {
        let p = Polygon::new(
            sq(0, 0, 10, 10),
            vec![Ring::from([(4, 4), (4, 6), (6, 6), (6, 4)])],
        );
        let t = trapezoids(&p).unwrap();
        // Left of the hole, below, above, right.
        assert_eq!(t.len(), 4);
        let a: f64 = t.iter().map(|t| t.area()).sum();
        assert_eq!(a, 96.0);
    }

    #[test]
    fn slanted() {
        let tri = Ring::from([(0, 0), (10, 1), (3, 7)]);
        let t = trapezoids(&tri).unwrap();
        assert_eq!(t.len(), 2);
        let a: f64 = t.iter().map(|t| t.area()).sum();
        assert!((a - 33.5).abs() < 1e-9);
    }
}
