//! Fracturing: polygons with holes to single hole-free outlines.
//!
//! Each hole is connected to the outline by a zero-width *cut-in*: a bridge segment walked
//! once in each direction. Holes are processed from left to right (by their
//! lexicographically smallest vertex `v`); a horizontal ray is cast from `v` towards `-x` and
//! the first outline edge it meets is found exactly. When that hit point is an integer point
//! the cut is horizontal (axis-aligned), splitting the hit edge if needed; otherwise the cut
//! goes to the outline vertex near the hit that is visible from `v` with the smallest angle
//! to the ray (as in the classic ear-clipping hole elimination). Holes touching the outline
//! at `v` are merged there without a bridge.
//!
//! All decisions use exact integer arithmetic; the result is deterministic.

use crate::error::{Result, check_point};
use crate::geom::{Point, Polygon, Rect, Ring};
use crate::predicates::{in_segment_interior, orient, segments_cross_properly};
use crate::query::ring_area2;
use core::cmp::Ordering;

/// Linked-list vertex of the outline under construction.
#[derive(Clone, Copy, Debug)]
struct Node {
    p: Point,
    next: u32,
    prev: u32,
    /// Ring index (0 = outer); used to know whether the node is part of the outline yet.
    ring: u32,
}

/// Uniform grid over node positions and edge bounding boxes (edges keyed by start node).
struct Grid {
    x0: i64,
    y0: i64,
    s: i64,
    nx: usize,
    ny: usize,
    edges: Vec<Vec<u32>>,
    nodes: Vec<Vec<u32>>,
}

impl Grid {
    fn col(&self, x: i64) -> usize {
        ((x - self.x0).div_euclid(self.s)).clamp(0, self.nx as i64 - 1) as usize
    }
    fn row(&self, y: i64) -> usize {
        ((y - self.y0).div_euclid(self.s)).clamp(0, self.ny as i64 - 1) as usize
    }
    fn add_edge(&mut self, id: u32, a: Point, b: Point) {
        let (c0, c1) = (self.col(a.x.min(b.x)), self.col(a.x.max(b.x)));
        let (r0, r1) = (self.row(a.y.min(b.y)), self.row(a.y.max(b.y)));
        for r in r0..=r1 {
            for c in c0..=c1 {
                self.edges[r * self.nx + c].push(id);
            }
        }
    }
    fn add_node(&mut self, id: u32, p: Point) {
        let c = self.row(p.y) * self.nx + self.col(p.x);
        self.nodes[c].push(id);
    }
}

struct Fracturer {
    nodes: Vec<Node>,
    merged: Vec<bool>,
    grid: Grid,
}

/// Rational `n / d` with `d > 0`.
#[derive(Clone, Copy, Debug)]
struct Q {
    n: i128,
    d: i128,
}

impl Q {
    fn cmp(self, o: Q) -> Ordering {
        (self.n * o.d).cmp(&(o.n * self.d))
    }
    fn int(v: i64) -> Q {
        Q { n: v as i128, d: 1 }
    }
}

impl Fracturer {
    fn edge(&self, i: u32) -> (Point, Point) {
        let n = &self.nodes[i as usize];
        (n.p, self.nodes[n.next as usize].p)
    }

    fn active(&self, i: u32) -> bool {
        self.merged[self.nodes[i as usize].ring as usize]
    }

    /// `true` when the direction from node `i` to `t` lies inside the outline's interior
    /// angle at `i` (interior on the left of the outline's edges).
    fn locally_inside(&self, i: u32, t: Point) -> bool {
        let n = &self.nodes[i as usize];
        let a = self.nodes[n.prev as usize].p;
        let p = n.p;
        let b = self.nodes[n.next as usize].p;
        if t == p {
            return true;
        }
        // Interior wedge: from direction (b - p) counter-clockwise to (a - p).
        if orient(a, p, b) >= 0 {
            // Convex (or straight) corner: strictly left of p->b and strictly left of a->p.
            orient(p, b, t) > 0 && orient(a, p, t) > 0
        } else {
            // Reflex corner: anything not in the exterior wedge.
            orient(p, b, t) > 0 || orient(a, p, t) > 0
        }
    }

    /// Whether outline node `k` can receive the bridge from hole node `vn`: the bridge
    /// direction must lie inside `k`'s interior angle, or, when both are at the same point,
    /// the hole's two edges must.
    fn accepts(&self, k: u32, vn: u32) -> bool {
        let v = &self.nodes[vn as usize];
        if self.nodes[k as usize].p == v.p {
            let a = self.nodes[v.prev as usize].p;
            let b = self.nodes[v.next as usize].p;
            return self.locally_inside(k, a) && self.locally_inside(k, b);
        }
        self.locally_inside(k, v.p)
    }

    /// `true` when segment `v-m` crosses no active outline edge and passes through no
    /// active vertex other than its endpoints.
    fn visible(&self, v: Point, m: Point) -> bool {
        let g = &self.grid;
        let (c0, c1) = (g.col(v.x.min(m.x)), g.col(v.x.max(m.x)));
        let (r0, r1) = (g.row(v.y.min(m.y)), g.row(v.y.max(m.y)));
        for r in r0..=r1 {
            for c in c0..=c1 {
                let cell = r * g.nx + c;
                for &e in &g.edges[cell] {
                    if !self.active(e) {
                        continue;
                    }
                    let (a, b) = self.edge(e);
                    if segments_cross_properly(v, m, a, b) {
                        return false;
                    }
                }
                for &k in &g.nodes[cell] {
                    if self.active(k) && in_segment_interior(v, m, self.nodes[k as usize].p) {
                        return false;
                    }
                }
            }
        }
        true
    }

    fn new_node(&mut self, p: Point, ring: u32) -> u32 {
        let id = self.nodes.len() as u32;
        self.nodes.push(Node {
            p,
            next: id,
            prev: id,
            ring,
        });
        self.grid.add_node(id, p);
        id
    }

    /// Splits edge `e` (from node `e` to its next) at `p`, returning the new node.
    fn split_edge(&mut self, e: u32, p: Point) -> u32 {
        let ring = self.nodes[e as usize].ring;
        let next = self.nodes[e as usize].next;
        let s = self.new_node(p, ring);
        self.nodes[s as usize].prev = e;
        self.nodes[s as usize].next = next;
        self.nodes[e as usize].next = s;
        self.nodes[next as usize].prev = s;
        let q = self.nodes[next as usize].p;
        self.grid.add_edge(s, p, q);
        s
    }

    /// Connects hole node `v` to outline node `m` with a zero-width bridge.
    fn bridge(&mut self, m: u32, v: u32) {
        let mp = self.nodes[m as usize].p;
        let vp = self.nodes[v as usize].p;
        let ring_m = self.nodes[m as usize].ring;
        let ring_v = self.nodes[v as usize].ring;
        if mp == vp {
            // Touching: splice the hole in at the shared point without a bridge:
            // a -> m -> h1 -> ... -> hk -> v -> b.
            let b = self.nodes[m as usize].next;
            let h1 = self.nodes[v as usize].next;
            let hk = self.nodes[v as usize].prev;
            self.nodes[m as usize].next = h1;
            self.nodes[h1 as usize].prev = m;
            self.nodes[hk as usize].next = v;
            self.nodes[v as usize].prev = hk;
            self.nodes[v as usize].next = b;
            self.nodes[b as usize].prev = v;
            // Edges keyed by m and v changed geometry: register them where they now lie.
            let (h1p, bp) = (self.nodes[h1 as usize].p, self.nodes[b as usize].p);
            self.grid.add_edge(m, mp, h1p);
            self.grid.add_edge(v, vp, bp);
            return;
        }
        let a = self.nodes[m as usize].prev;
        let hk = self.nodes[v as usize].prev;
        // a -> m' -> v -> ... -> hk -> v' -> m -> b
        let m2 = self.new_node(mp, ring_m);
        let v2 = self.new_node(vp, ring_v);
        self.nodes[a as usize].next = m2;
        self.nodes[m2 as usize].prev = a;
        self.nodes[m2 as usize].next = v;
        self.nodes[v as usize].prev = m2;
        self.nodes[hk as usize].next = v2;
        self.nodes[v2 as usize].prev = hk;
        self.nodes[v2 as usize].next = m;
        self.nodes[m as usize].prev = v2;
        self.grid.add_edge(m2, mp, vp);
        self.grid.add_edge(v2, vp, mp);
    }

    /// Finds the outline node to bridge hole vertex `v` to.
    fn find_bridge(&mut self, vn: u32) -> Option<u32> {
        let v = self.nodes[vn as usize].p;
        // 1. Ray towards -x: nearest active edge going downwards (interior on its east side)
        //    with p.y >= v.y >= q.y, p.y > q.y.
        let mut best: Option<(Q, u32)> = None;
        let row = self.grid.row(v.y);
        let mut col = self.grid.col(v.x) as i64;
        while col >= 0 {
            let cell_min_x = self.grid.x0 + col * self.grid.s;
            if let Some((bx, _)) = best {
                // Every hit in this cell is at x <= cell max x; stop when that is left of best.
                if Q::int(cell_min_x + self.grid.s).cmp(bx) == Ordering::Less {
                    break;
                }
            }
            let cell = row * self.grid.nx + col as usize;
            for k in 0..self.grid.edges[cell].len() {
                let e = self.grid.edges[cell][k];
                if !self.active(e) {
                    continue;
                }
                let (p, q) = self.edge(e);
                if !(p.y >= v.y && v.y >= q.y && p.y > q.y) {
                    continue;
                }
                let dy = (p.y - q.y) as i128;
                let x = Q {
                    n: p.x as i128 * dy + (p.y - v.y) as i128 * (q.x - p.x) as i128,
                    d: dy,
                };
                if x.cmp(Q::int(v.x)) == Ordering::Greater {
                    continue;
                }
                let better = match best {
                    None => true,
                    Some((bx, be)) => match x.cmp(bx) {
                        Ordering::Greater => true,
                        // Same hit point (shared vertex): keep the smaller node id.
                        Ordering::Equal => e < be,
                        Ordering::Less => false,
                    },
                };
                if better {
                    best = Some((x, e));
                }
            }
            col -= 1;
        }
        let (hx, he) = best?;
        let (p, q) = self.edge(he);
        // 2. Integer hit point: horizontal cut (or direct merge when touching).
        if hx.n % hx.d == 0 {
            let h = Point::new((hx.n / hx.d) as i64, v.y);
            let target = if h == p {
                Some(he)
            } else if h == q {
                Some(self.nodes[he as usize].next)
            } else {
                None
            };
            let t = match target {
                Some(t) => {
                    // Several nodes may sit at h (earlier bridges, touching rings): pick one
                    // whose interior wedge contains the direction towards v.
                    self.node_at(h, vn).unwrap_or(t)
                }
                None => self.split_edge(he, h),
            };
            return Some(t);
        }
        // 3. Rational hit: the endpoint of the hit edge with the smaller x, or a vertex in
        //    the triangle (v, hit, m) with the smallest angle to the ray.
        let mut m = if p.x < q.x {
            he
        } else {
            self.nodes[he as usize].next
        };
        let mp = self.nodes[m as usize].p;
        // Triangle bbox.
        let hx_floor = hx.n.div_euclid(hx.d) as i64;
        let tri = Rect::new(
            Point::new(mp.x.min(hx_floor), mp.y.min(v.y)),
            Point::new(v.x, mp.y.max(v.y)),
        );
        let mut best_tan: Option<(i128, i128)> = None; // |dy| / dx
        let g = &self.grid;
        let (c0, c1) = (g.col(tri.min.x), g.col(tri.max.x));
        let (r0, r1) = (g.row(tri.min.y), g.row(tri.max.y));
        let mut cands: Vec<u32> = Vec::new();
        for r in r0..=r1 {
            for c in c0..=c1 {
                cands.extend(
                    g.nodes[r * g.nx + c]
                        .iter()
                        .copied()
                        .filter(|&k| self.active(k)),
                );
            }
        }
        cands.sort_unstable();
        for k in cands {
            let pp = self.nodes[k as usize].p;
            if !(pp.x <= v.x && pp.x >= mp.x && pp.x != v.x) {
                continue;
            }
            // Inside the closed triangle (v, hit, m)? Use orientation tests; the hit point is
            // rational, so scale: hit = (hx.n / hx.d, v.y).
            if !in_triangle_hit(v, hx, mp, pp) {
                continue;
            }
            if !self.locally_inside(k, v) {
                continue;
            }
            let dy = (v.y - pp.y).unsigned_abs() as i128;
            let dx = (v.x - pp.x) as i128;
            let better = match best_tan {
                None => true,
                Some((by, bx)) => {
                    let o = (dy * bx).cmp(&(by * dx));
                    o == Ordering::Less
                        || (o == Ordering::Equal && pp.x > self.nodes[m as usize].p.x)
                }
            };
            if better {
                best_tan = Some((dy, dx));
                m = k;
            }
        }
        Some(m)
    }

    /// An active node at `h` that [`accepts`](Self::accepts) the bridge from `vn`.
    fn node_at(&self, h: Point, vn: u32) -> Option<u32> {
        let g = &self.grid;
        let cell = g.row(h.y) * g.nx + g.col(h.x);
        let mut ks: Vec<u32> = g.nodes[cell]
            .iter()
            .copied()
            .filter(|&k| self.active(k) && self.nodes[k as usize].p == h)
            .collect();
        ks.sort_unstable();
        ks.into_iter().find(|&k| self.accepts(k, vn))
    }

    /// Fallback: any visible, locally-inside active node, nearest first.
    fn brute_bridge(&self, vn: u32) -> Option<u32> {
        let v = self.nodes[vn as usize].p;
        let mut cands: Vec<(i128, u32)> = (0..self.nodes.len() as u32)
            .filter(|&k| self.active(k))
            .map(|k| (crate::predicates::dist2(v, self.nodes[k as usize].p), k))
            .collect();
        cands.sort_unstable();
        cands
            .into_iter()
            .map(|c| c.1)
            .find(|&k| self.accepts(k, vn) && self.visible(v, self.nodes[k as usize].p))
    }
}

/// `true` when `p` lies in the closed triangle `(v, hit, m)` with `hit = (hx, v.y)` rational.
fn in_triangle_hit(v: Point, hx: Q, m: Point, p: Point) -> bool {
    // Scale everything by hx.d so the hit becomes integer (values stay within i128).
    let d = hx.d;
    let s = |q: Point| (q.x as i128 * d, q.y as i128 * d);
    let (vx, vy) = s(v);
    let (mx, my) = s(m);
    let (px, py) = s(p);
    let (hxx, hyy) = (hx.n, v.y as i128 * d);
    let or = |ax: i128, ay: i128, bx: i128, by: i128, cx: i128, cy: i128| -> i128 {
        // Values up to 2^42 * 2^42 scaled products: compare signs via f64-safe split.
        let l = (bx - ax) as f64 * (cy - ay) as f64;
        let r = (by - ay) as f64 * (cx - ax) as f64;
        let diff = l - r;
        let mag = l.abs() + r.abs();
        if diff.abs() > mag * 1e-12 {
            return if diff > 0.0 { 1 } else { -1 };
        }
        // Exact fallback with reduced values (divide out the common scale where possible).
        exact_orient_sign(ax, ay, bx, by, cx, cy)
    };
    let o1 = or(vx, vy, hxx, hyy, px, py);
    let o2 = or(hxx, hyy, mx, my, px, py);
    let o3 = or(mx, my, vx, vy, px, py);
    (o1 >= 0 && o2 >= 0 && o3 >= 0) || (o1 <= 0 && o2 <= 0 && o3 <= 0)
}

/// Exact sign of `(b - a) x (c - a)` for large `i128` inputs, via 256-bit products.
fn exact_orient_sign(ax: i128, ay: i128, bx: i128, by: i128, cx: i128, cy: i128) -> i128 {
    let mul = |a: i128, b: i128| -> (bool, crate::wide::U384) {
        let neg = (a < 0) != (b < 0) && a != 0 && b != 0;
        (
            neg,
            crate::wide::U384::mul_u128(a.unsigned_abs(), b.unsigned_abs()),
        )
    };
    let (n1, l) = mul(bx - ax, cy - ay);
    let (n2, r) = mul(by - ay, cx - ax);
    // sign(l_signed - r_signed)
    match (n1, n2) {
        (false, true) => {
            if l == crate::wide::U384::default() && r == crate::wide::U384::default() {
                0
            } else {
                1
            }
        }
        (true, false) => {
            if l == crate::wide::U384::default() && r == crate::wide::U384::default() {
                0
            } else {
                -1
            }
        }
        (false, false) => match l.cmp(&r) {
            Ordering::Greater => 1,
            Ordering::Less => -1,
            Ordering::Equal => 0,
        },
        (true, true) => match l.cmp(&r) {
            Ordering::Greater => -1,
            Ordering::Less => 1,
            Ordering::Equal => 0,
        },
    }
}

/// Fractures a polygon into a single outline: its holes are joined to the outer ring by
/// zero-width cut-ins, as required for Gerber regions.
///
/// The result is counter-clockwise and *weakly* simple: every cut is walked once in each
/// direction, so cut endpoints appear twice. Its signed area equals the polygon's area, and
/// the region it encloses (non-zero rule) is exactly the polygon. Cut placement is
/// deterministic; cuts are horizontal whenever the horizontal ray from a hole's
/// lexicographically smallest vertex meets the outline at an integer point. Consecutive
/// duplicate vertices (from holes touching the outline) are removed.
///
/// Rings are normalized first (outer counter-clockwise, holes clockwise); the input should
/// be a valid polygon (see [`validate`](crate::validate)) — for invalid input the result is
/// unspecified but no panic occurs.
///
/// ```
/// use polyclip::{fracture, Polygon, Ring};
/// let p = Polygon::new(
///     Ring::from([(0, 0), (10, 0), (10, 10), (0, 10)]),
///     vec![Ring::from([(4, 4), (4, 6), (6, 6), (6, 4)])],
/// );
/// let f = fracture(&p).unwrap();
/// assert_eq!(f.signed_area2(), p.signed_area2());
/// assert!(f.contains(&polyclip::Point::new(0, 4))); // horizontal cut from (4, 4) to (0, 4)
/// ```
pub fn fracture(poly: &Polygon) -> Result<Ring> {
    for r in poly.rings() {
        for &p in r.iter() {
            check_point(p)?;
        }
    }
    let clean = |r: &Ring, ccw: bool| -> Vec<Point> {
        let mut v: Vec<Point> = Vec::with_capacity(r.len());
        for &p in r.iter() {
            if v.last() != Some(&p) {
                v.push(p);
            }
        }
        while v.len() > 1 && v.first() == v.last() {
            v.pop();
        }
        if (ring_area2(&v) > 0) != ccw {
            v.reverse();
        }
        v
    };
    let outer = clean(&poly.outer, true);
    if outer.len() < 3 || ring_area2(&outer) == 0 {
        return Ok(Ring(outer));
    }
    let mut rings: Vec<Vec<Point>> = vec![outer];
    for h in &poly.holes {
        let h = clean(h, false);
        if h.len() >= 3 && ring_area2(&h) != 0 {
            rings.push(h);
        }
    }
    if rings.len() == 1 {
        return Ok(canonical_start(rings.pop().unwrap()));
    }
    // Grid sized from the total extent and vertex count.
    let all = Rect::of_points(rings.iter().flatten()).unwrap();
    let nv: usize = rings.iter().map(|r| r.len()).sum();
    let w = (all.width() + 1) as f64;
    let h = (all.height() + 1) as f64;
    let s = (libm::ceil(libm::sqrt(w * h / (nv as f64 / 2.0 + 1.0))) as i64).max(1);
    let nx = (all.width() / s + 1) as usize;
    let ny = (all.height() / s + 1) as usize;
    let grid = Grid {
        x0: all.min.x,
        y0: all.min.y,
        s,
        nx,
        ny,
        edges: vec![Vec::new(); nx * ny],
        nodes: vec![Vec::new(); nx * ny],
    };
    let mut f = Fracturer {
        nodes: Vec::with_capacity(nv * 2),
        merged: vec![false; rings.len()],
        grid,
    };
    f.merged[0] = true;
    let mut first_node: Vec<u32> = Vec::with_capacity(rings.len());
    for (ri, r) in rings.iter().enumerate() {
        let base = f.nodes.len() as u32;
        let n = r.len() as u32;
        first_node.push(base);
        for (k, &p) in r.iter().enumerate() {
            let k = k as u32;
            f.nodes.push(Node {
                p,
                next: base + (k + 1) % n,
                prev: base + (k + n - 1) % n,
                ring: ri as u32,
            });
        }
        for k in 0..n {
            let id = base + k;
            f.grid.add_node(id, r[k as usize]);
            f.grid
                .add_edge(id, r[k as usize], r[((k + 1) % n) as usize]);
        }
    }
    // Holes by lexicographically smallest vertex.
    let mut order: Vec<(Point, u32, usize)> = (1..rings.len())
        .map(|ri| {
            let (k, p) = rings[ri]
                .iter()
                .enumerate()
                .min_by_key(|(_, p)| **p)
                .unwrap();
            (*p, first_node[ri] + k as u32, ri)
        })
        .collect();
    order.sort_unstable();
    for (_, vn, ri) in order {
        let v = f.nodes[vn as usize].p;
        let m = f
            .find_bridge(vn)
            .filter(|&m| f.accepts(m, vn) && f.visible(v, f.nodes[m as usize].p))
            .or_else(|| f.brute_bridge(vn));
        let Some(m) = m else {
            // Unbridgeable (invalid input): leave the hole out.
            continue;
        };
        f.bridge(m, vn);
        f.merged[ri] = true;
    }
    // Walk the outline from the outer ring's smallest vertex.
    let start = {
        let (k, _) = rings[0]
            .iter()
            .enumerate()
            .min_by_key(|(_, p)| **p)
            .unwrap();
        first_node[0] + k as u32
    };
    let mut out: Vec<Point> = Vec::with_capacity(f.nodes.len());
    let mut k = start;
    loop {
        let p = f.nodes[k as usize].p;
        if out.last() != Some(&p) {
            out.push(p);
        }
        k = f.nodes[k as usize].next;
        if k == start || out.len() > f.nodes.len() {
            break;
        }
    }
    while out.len() > 1 && out.first() == out.last() {
        out.pop();
    }
    Ok(Ring(out))
}

fn canonical_start(mut v: Vec<Point>) -> Ring {
    if let Some((i, _)) = v.iter().enumerate().min_by_key(|(_, p)| **p) {
        v.rotate_left(i);
    }
    Ring(v)
}

/// Fractures every polygon of a set (see [`fracture`]).
pub fn fracture_set(polys: &[Polygon]) -> Result<Vec<Ring>> {
    polys.iter().map(fracture).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{FillRule, union_all};

    fn sq(x0: i64, y0: i64, x1: i64, y1: i64) -> Ring {
        Ring::from([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])
    }

    /// The fractured ring must enclose exactly the polygon and never cross itself.
    fn check(p: &Polygon) {
        let f = fracture(p).unwrap();
        let orig = union_all(p, FillRule::EvenOdd).unwrap();
        let a: i128 = orig.iter().map(|q| q.signed_area2()).sum();
        assert_eq!(f.signed_area2(), a, "{f:?}");
        let n = f.len();
        for i in 0..n {
            for j in i + 1..n {
                let (a, b) = (f[i], f[(i + 1) % n]);
                let (c, d) = (f[j], f[(j + 1) % n]);
                assert!(
                    !segments_cross_properly(a, b, c, d),
                    "crossing {a:?}-{b:?} {c:?}-{d:?}"
                );
            }
        }
        let back = union_all(&f, FillRule::NonZero).unwrap();
        assert_eq!(back, orig);
    }

    #[test]
    fn single_hole_horizontal_cut() {
        let p = Polygon::new(sq(0, 0, 10, 10), vec![sq(4, 4, 6, 6)]);
        let f = fracture(&p).unwrap();
        assert_eq!(
            f,
            Ring::from([
                (0, 0),
                (10, 0),
                (10, 10),
                (0, 10),
                (0, 4),
                (4, 4),
                (4, 6),
                (6, 6),
                (6, 4),
                (4, 4),
                (0, 4)
            ])
        );
        check(&p);
    }

    #[test]
    fn many_holes() {
        let mut holes = Vec::new();
        for i in 0..5 {
            for j in 0..5 {
                holes.push(sq(2 + 4 * i, 2 + 4 * j, 4 + 4 * i, 4 + 4 * j));
            }
        }
        check(&Polygon::new(sq(0, 0, 22, 22), holes));
    }

    #[test]
    fn rational_hit_and_touching() {
        // Slanted left side: ray from (5, 5) hits at x = 2.5.
        let outer = Ring::from([(0, 0), (20, 0), (20, 20), (5, 20)]);
        check(&Polygon::new(outer.clone(), vec![sq(6, 5, 8, 7)]));
        // Hole touching the outer ring at its smallest vertex.
        let p = union_all(&vec![sq(0, 0, 10, 10)], FillRule::NonZero).unwrap();
        let d = crate::boolean(
            crate::Op::Difference,
            &p,
            &Ring::from([(0, 5), (5, 3), (5, 7)]),
            FillRule::NonZero,
        )
        .unwrap();
        check(&d[0]);
        // Hole touching another hole.
        let d = crate::boolean(
            crate::Op::Difference,
            &sq(0, 0, 20, 20),
            &vec![sq(2, 2, 6, 6), sq(6, 6, 10, 10)],
            FillRule::NonZero,
        )
        .unwrap();
        assert_eq!(d[0].holes.len(), 2);
        check(&d[0]);
    }
}
