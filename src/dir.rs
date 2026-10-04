//! Exact integer projection directions for sweep-and-prune.

use crate::geom::Point;

/// A projection direction for sweep-and-prune: the exact integer projection
/// `p.x * nx + p.y * ny`. Candidates are the axes, the diagonals and the normal of the
/// longest segment of the leaf, so dense parallel segments at any angle separate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Dir {
    pub nx: i64,
    pub ny: i64,
}

impl Dir {
    pub const X: Dir = Dir { nx: 1, ny: 0 };
    pub(crate) const FIXED: [Dir; 4] = [
        Dir { nx: 1, ny: 0 },
        Dir { nx: 0, ny: 1 },
        Dir { nx: 1, ny: 1 },
        Dir { nx: 1, ny: -1 },
    ];

    #[inline]
    pub fn proj(self, p: Point) -> i128 {
        p.x as i128 * self.nx as i128 + p.y as i128 * self.ny as i128
    }

    /// Projected interval of a segment.
    #[inline]
    pub fn range(self, s: &(Point, Point)) -> (i128, i128) {
        let (a, b) = (self.proj(s.0), self.proj(s.1));
        (a.min(b), a.max(b))
    }

    /// Projected margin covering distance 1 (with slack): `ceil(|n|) + 1`.
    #[inline]
    pub fn margin(self) -> i128 {
        libm::ceil(libm::hypot(self.nx as f64, self.ny as f64)) as i128 + 1
    }

    /// The direction along which the segments' projections are thinnest in total
    /// (projected widths are compared after dividing by `|n|`).
    pub fn best(segs: &[(Point, Point)], items: &[u32]) -> Dir {
        Self::best_of(items.iter().map(|&i| segs[i as usize]))
    }

    /// Like [`best`](Self::best) over any segment iterator (iterated twice).
    pub fn best_of<I: Iterator<Item = (Point, Point)> + Clone>(segs: I) -> Dir {
        // Normal of the longest segment (reduced), as a data-driven candidate.
        let longest = segs
            .clone()
            .max_by_key(|s| crate::predicates::dist2(s.0, s.1))
            .map(|s| {
                let (dx, dy) = (s.1.x - s.0.x, s.1.y - s.0.y);
                let g = gcd(dx.unsigned_abs(), dy.unsigned_abs()).max(1) as i64;
                Dir {
                    nx: -dy / g,
                    ny: dx / g,
                }
            });
        let mut best = (f64::INFINITY, Dir::X);
        for d in Dir::FIXED.into_iter().chain(longest) {
            if d.nx == 0 && d.ny == 0 {
                continue;
            }
            let w: i128 = segs
                .clone()
                .map(|s| {
                    let (lo, hi) = d.range(&s);
                    hi - lo
                })
                .sum();
            let w = w as f64 / libm::hypot(d.nx as f64, d.ny as f64);
            if w < best.0 {
                best = (w, d);
            }
        }
        best.1
    }
}

fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

/// Direction for searching pairs between two sets: the one along which the sets' projected
/// extents overlap least (relative to their union); among separating directions the one
/// with the widest gap, otherwise the one with the thinnest segments.
pub(crate) fn separating(
    a: impl Iterator<Item = (Point, Point)> + Clone,
    b: impl Iterator<Item = (Point, Point)> + Clone,
) -> Dir {
    let thin = Dir::best_of(a.clone().chain(b.clone()));
    let mut best = ((f64::INFINITY, f64::INFINITY, f64::INFINITY), Dir::X);
    for d in Dir::FIXED.into_iter().chain([thin]) {
        let ext = |it: &mut dyn Iterator<Item = (Point, Point)>| -> (i128, i128, i128) {
            let (mut lo, mut hi, mut w) = (i128::MAX, i128::MIN, 0i128);
            for s in it {
                let (l, h) = d.range(&s);
                lo = lo.min(l);
                hi = hi.max(h);
                w += h - l;
            }
            (lo, hi, w)
        };
        let (alo, ahi, aw) = ext(&mut a.clone());
        let (blo, bhi, bw) = ext(&mut b.clone());
        let inner = (ahi.min(bhi) - alo.max(blo)) as f64;
        let union = (ahi.max(bhi) - alo.min(blo)).max(1) as f64;
        let norm = libm::hypot(d.nx as f64, d.ny as f64);
        // (overlap ratio, minus the normalized gap, normalized total width): smaller is better.
        let score = (
            inner.max(0.0) / union,
            (inner.min(0.0)) / norm,
            (aw + bw) as f64 / norm,
        );
        if score < best.0 {
            best = (score, d);
        }
    }
    best.1
}
