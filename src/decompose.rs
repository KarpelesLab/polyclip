//! Trapezoidal (vertical) decomposition.

use crate::boolean::{FillRule, RingSource, union_all};
use crate::error::Result;
use crate::geom::{Point, PointF};
use crate::sweep::{cmp_sweep_edges, sweep_events};

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

    /// Height at `x` (top minus bottom), from an exact rational evaluated once in `f64`.
    fn height(&self, x: i64) -> f64 {
        let ((a, b), (c, d)) = (self.top, self.bottom);
        let (dxt, dyt) = ((b.x - a.x) as i128, (b.y - a.y) as i128);
        let (dxb, dyb) = ((d.x - c.x) as i128, (d.y - c.y) as i128);
        // Relative to the bottom edge's start to keep the numbers small.
        let (ox, oy) = (c.x as i128, c.y as i128);
        let top = ((a.y as i128 - oy) * dxt + (x as i128 - a.x as i128) * dyt) * dxb;
        let bot = (x as i128 - ox) * dyb * dxt;
        (top - bot) as f64 / (dxt * dxb) as f64
    }

    /// Area (as `f64`, computed from exact heights).
    pub fn area(&self) -> f64 {
        let w = (self.x1 - self.x0) as f64;
        w * (self.height(self.x0) + self.height(self.x1)) / 2.0
    }
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
    // Boundary edges (canonical rings: interior on the left), as (lo, hi) with the
    // interior above when the ring runs from lo to hi.
    let mut edges: Vec<(Point, Point, bool)> = Vec::new();
    for p in &polys {
        for r in p.rings() {
            for (a, b) in r.edges() {
                edges.push(if a < b { (a, b, true) } else { (b, a, false) });
            }
        }
    }
    edges.sort_unstable_by(|a, b| cmp_sweep_edges((a.0, a.1), (b.0, b.1)));
    let segs: Vec<(Point, Point)> = edges.iter().map(|e| (e.0, e.1)).collect();
    // For every edge with the interior above it: the trapezoid open in that gap
    // (start x, top edge).
    let mut open: Vec<Option<(i64, u32)>> = vec![None; edges.len()];
    let mut out: Vec<Trapezoid> = Vec::new();
    let emit = |bottom: u32, x0: i64, top: u32, x1: i64, out: &mut Vec<Trapezoid>| {
        // Zero-width pieces (between walls at the same x, or along vertical edges) vanish.
        if x1 > x0 {
            let (b, t) = (&edges[bottom as usize], &edges[top as usize]);
            out.push(Trapezoid {
                x0,
                x1,
                bottom: (b.0, b.1),
                top: (t.0, t.1),
            });
        }
    };
    sweep_events(&segs, |below, above, ending, starting| {
        let v = if starting.is_empty() {
            segs[ending[0] as usize].1
        } else {
            segs[starting.start as usize].0
        };
        let x = v.x;
        // A wall through `v` closes the trapezoid of every gap touching the vertex.
        if let Some(b) = below
            && let Some((x0, t)) = open[b as usize].take()
        {
            emit(b, x0, t, x, &mut out);
        }
        for &e in ending {
            if let Some((x0, t)) = open[e as usize].take() {
                emit(e, x0, t, x, &mut out);
            }
        }
        // New gaps, bottom to top: below's gap up to the first new edge (or `above`), the
        // gaps between new edges, and the gap above the last new edge.
        let mut bottoms: Vec<u32> = below.into_iter().collect();
        bottoms.extend(starting.clone());
        for (k, &bot) in bottoms.iter().enumerate() {
            let top = bottoms.get(k + 1).copied().or(above);
            if let Some(top) = top
                && edges[bot as usize].2
            {
                open[bot as usize] = Some((x, top));
            }
        }
    });
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
