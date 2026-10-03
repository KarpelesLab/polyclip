//! Shared input decoding for the polyclip fuzz targets.
//!
//! Every target decodes its raw libFuzzer bytes through [`Gen`], which turns them into
//! geometry biased towards the cases that break polygon code:
//!
//! * a per-input (and occasionally per-ring) coordinate **scale**: a 4×4 grid (heavy
//!   degeneracy: almost every edge pair is collinear, touching or coincident), small,
//!   medium, the full `±MAX_COORD` range, and values hugging the range boundary;
//! * rarely (only when the input enables it) coordinates **just beyond** `MAX_COORD` or at
//!   `i64::MIN`/`i64::MAX`, which every operation must reject with `Err`, never panic;
//! * **duplicate** points, **collinear** continuations of the previous edge, points reused
//!   from a shared pool (coincident vertices across rings), axis-aligned steps;
//! * CAD-like rings: rectangles and approximated circles mixed with arbitrary rings.
//!
//! Decoding never fails: when the bytes run out, `arbitrary` returns the low end of each
//! range, so short inputs decode to small, well-formed geometry.

use arbitrary::Unstructured;
use polyclip::{
    ArcTol, Circle, EndCap, FillRule, Join, MAX_COORD, Op, Path, Point, Polygon, Ring, Side,
};

/// Coordinate scale.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Scale {
    /// Coordinates in `0..=3`: maximal degeneracy.
    Tiny,
    /// Coordinates in `-16..=16`.
    Small,
    /// Coordinates in `±2^20`.
    Medium,
    /// Coordinates anywhere in `±MAX_COORD`.
    Full,
    /// Coordinates at or next to `±MAX_COORD` (and 0).
    Edge,
}

/// Byte-driven geometry generator.
pub struct Gen<'a> {
    /// Underlying byte source.
    pub u: Unstructured<'a>,
    scale: Scale,
    /// Whether out-of-range coordinates may be produced.
    allow_oor: bool,
    /// Set once an out-of-range coordinate has been produced.
    pub produced_oor: bool,
    pool: Vec<Point>,
}

const M: i64 = MAX_COORD;

impl<'a> Gen<'a> {
    /// Starts decoding `data`; the first bytes select the scale and the out-of-range mode.
    pub fn new(data: &'a [u8]) -> Self {
        let mut u = Unstructured::new(data);
        let scale = pick_scale(&mut u);
        // One input in 16 may contain out-of-range coordinates.
        let allow_oor = u.int_in_range(0u8..=15).unwrap_or(0) == 15;
        Gen {
            u,
            scale,
            allow_oor,
            produced_oor: false,
            pool: Vec::new(),
        }
    }

    /// Current scale.
    pub fn scale(&self) -> Scale {
        self.scale
    }

    /// A byte (0 when exhausted).
    pub fn byte(&mut self) -> u8 {
        self.u.int_in_range(0u8..=255).unwrap_or(0)
    }

    /// An integer in `lo..=hi` (`lo` when exhausted).
    pub fn int(&mut self, lo: i64, hi: i64) -> i64 {
        self.u.int_in_range(lo..=hi).unwrap_or(lo)
    }

    /// A boolean.
    pub fn bool(&mut self) -> bool {
        self.byte() & 1 == 1
    }

    /// Typical magnitude of coordinates at the current scale (used to size deltas,
    /// radii and tolerances).
    pub fn extent(&self) -> i64 {
        match self.scale {
            Scale::Tiny => 4,
            Scale::Small => 16,
            Scale::Medium => 1 << 20,
            Scale::Full | Scale::Edge => M,
        }
    }

    /// One coordinate at the current scale.
    pub fn coord(&mut self) -> i64 {
        if self.allow_oor && self.int(0, 31) == 0 {
            self.produced_oor = true;
            let v = *self
                .u
                .choose(&[M + 1, -M - 1, M + 2, 2 * M, i64::MAX, i64::MIN, -i64::MAX])
                .unwrap_or(&(M + 1));
            return v;
        }
        match self.scale {
            Scale::Tiny => self.int(0, 3),
            Scale::Small => self.int(-16, 16),
            Scale::Medium => self.int(-(1 << 20), 1 << 20),
            Scale::Full => self.int(-M, M),
            Scale::Edge => {
                let base = *self.u.choose(&[M, -M, 0]).unwrap_or(&M);
                let j = self.int(0, 3);
                // Move towards the inside of the range.
                if base > 0 {
                    base - j
                } else if base < 0 {
                    base + j
                } else {
                    j - 1
                }
            }
        }
    }

    /// A fresh point.
    pub fn point(&mut self) -> Point {
        Point::new(self.coord(), self.coord())
    }

    fn note(&mut self, p: Point) -> Point {
        if !p.in_range() {
            self.produced_oor = true;
        }
        if self.pool.len() < 64 {
            self.pool.push(p);
        }
        p
    }

    /// Next vertex of a polyline whose vertices so far are `prev`.
    fn next_point(&mut self, prev: &[Point]) -> Point {
        let k = self.int(0, 9);
        let p = match (k, prev) {
            // Duplicate of the previous point.
            (0, [.., l]) => *l,
            // Collinear with the previous edge: a + t (b - a), t in -2..=3.
            (1 | 2, [.., a, b]) => {
                let t = self.int(-2, 3) as i128;
                let x = a.x as i128 + t * (b.x as i128 - a.x as i128);
                let y = a.y as i128 + t * (b.y as i128 - a.y as i128);
                let c = |v: i128| v.clamp(-(M as i128), M as i128) as i64;
                Point::new(c(x), c(y))
            }
            // A point used before (coincident vertices across rings / self-touch).
            (3, _) if !self.pool.is_empty() => {
                let i = self.int(0, self.pool.len() as i64 - 1) as usize;
                self.pool[i]
            }
            // Axis-aligned step.
            (4, [.., l]) => {
                let c = self.coord();
                if self.bool() {
                    Point::new(c, l.y)
                } else {
                    Point::new(l.x, c)
                }
            }
            _ => self.point(),
        };
        self.note(p)
    }

    /// An arbitrary ring (possibly degenerate, self-intersecting, any orientation).
    pub fn raw_ring(&mut self, max_pts: usize) -> Ring {
        let n = self.int(0, max_pts as i64) as usize;
        let mut v: Vec<Point> = Vec::with_capacity(n);
        for _ in 0..n {
            let p = self.next_point(&v);
            v.push(p);
        }
        Ring(v)
    }

    /// An arbitrary open path.
    pub fn path(&mut self, max_pts: usize) -> Path {
        Path(self.raw_ring(max_pts).0)
    }

    /// An axis-aligned rectangle, either orientation.
    pub fn rect(&mut self) -> Ring {
        let a = self.point();
        let b = self.point();
        self.note(a);
        self.note(b);
        let mut r = Ring::from([(a.x, a.y), (b.x, a.y), (b.x, b.y), (a.x, b.y)]);
        if self.bool() {
            r.reverse_orientation();
        }
        r
    }

    /// A circle approximated by polyclip itself (falls back to a rectangle when the
    /// circle is out of range).
    pub fn circle_ring(&mut self) -> Ring {
        let c = self.point();
        let ext = self.extent().max(4);
        let r = self.int(1, ext / 2 + 1);
        let tol = self.arc_tol_for(r);
        match Circle::new(c, r).to_ring(tol) {
            Ok(mut ring) => {
                if self.bool() {
                    ring.reverse_orientation();
                }
                ring
            }
            Err(_) => self.rect(),
        }
    }

    /// A ring of any kind; occasionally switches the coordinate scale first.
    pub fn ring(&mut self, max_pts: usize) -> Ring {
        if self.int(0, 7) == 0 {
            self.scale = pick_scale(&mut self.u);
        }
        match self.int(0, 7) {
            0 => self.rect(),
            1 => self.circle_ring(),
            _ => self.raw_ring(max_pts),
        }
    }

    /// Up to `max_rings` rings.
    pub fn rings(&mut self, max_rings: usize, max_pts: usize) -> Vec<Ring> {
        let n = self.int(0, max_rings as i64) as usize;
        (0..n).map(|_| self.ring(max_pts)).collect()
    }

    /// Up to `max_paths` paths.
    pub fn paths(&mut self, max_paths: usize, max_pts: usize) -> Vec<Path> {
        let n = self.int(0, max_paths as i64) as usize;
        (0..n).map(|_| self.path(max_pts)).collect()
    }

    /// Polygons built from raw rings (outer + holes, not normalized).
    pub fn raw_polygons(&mut self, max_polys: usize) -> Vec<Polygon> {
        let n = self.int(0, max_polys as i64) as usize;
        (0..n)
            .map(|_| {
                let outer = self.ring(10);
                let nh = self.int(0, 3) as usize;
                let holes = (0..nh).map(|_| self.ring(8)).collect();
                Polygon { outer, holes }
            })
            .collect()
    }

    /// A fill rule.
    pub fn fill_rule(&mut self) -> FillRule {
        *self
            .u
            .choose(&[
                FillRule::EvenOdd,
                FillRule::NonZero,
                FillRule::Positive,
                FillRule::Negative,
            ])
            .unwrap_or(&FillRule::NonZero)
    }

    /// A boolean operation.
    pub fn op(&mut self) -> Op {
        *self
            .u
            .choose(&[Op::Union, Op::Intersection, Op::Difference, Op::Xor])
            .unwrap_or(&Op::Union)
    }

    /// A join type (miter limits include degenerate values: < 1, NaN, huge).
    pub fn join(&mut self) -> Join {
        match self.int(0, 3) {
            0 => Join::Round,
            1 => Join::Bevel,
            2 => Join::Square,
            _ => {
                let limit = *self
                    .u
                    .choose(&[0.0, 0.5, 1.0, 1.5, 2.0, 4.0, 100.0, f64::NAN, f64::INFINITY])
                    .unwrap_or(&2.0);
                Join::Miter { limit }
            }
        }
    }

    /// An end cap.
    pub fn cap(&mut self) -> EndCap {
        *self
            .u
            .choose(&[EndCap::Round, EndCap::Square, EndCap::Butt, EndCap::Joined])
            .unwrap_or(&EndCap::Round)
    }

    /// An arc side.
    pub fn side(&mut self) -> Side {
        *self
            .u
            .choose(&[Side::Outside, Side::Inside, Side::Nearest])
            .unwrap_or(&Side::Outside)
    }

    /// An arc tolerance suited to radius `r`: never so fine that a full circle needs more
    /// than ~ a few hundred vertices (the fuzzer would otherwise spend its time generating
    /// millions of arc points), occasionally invalid (`< 1`) to exercise the error path.
    pub fn arc_tol_for(&mut self, r: i64) -> ArcTol {
        let side = self.side();
        if self.int(0, 63) == 0 {
            return ArcTol::new(self.int(-2, 0), side);
        }
        let floor = (r.unsigned_abs() >> 14).max(1) as i64;
        let t = self.int(1, 64).max(floor);
        ArcTol::new(t, side)
    }

    /// An offset distance at the current scale (either sign).
    pub fn delta(&mut self) -> i64 {
        let ext = self.extent();
        match self.int(0, 7) {
            0 => 0,
            1 => self.int(-3, 3),
            2 => *self
                .u
                .choose(&[2 * M, -2 * M, 2 * M + 1, i64::MAX, i64::MIN])
                .unwrap_or(&0),
            _ => self.int(-ext, ext),
        }
    }
}

fn pick_scale(u: &mut Unstructured<'_>) -> Scale {
    *u.choose(&[
        Scale::Tiny,
        Scale::Tiny,
        Scale::Small,
        Scale::Small,
        Scale::Medium,
        Scale::Full,
        Scale::Edge,
    ])
    .unwrap_or(&Scale::Tiny)
}

/// `true` when every point is within `±MAX_COORD`.
pub fn all_in_range<'p>(pts: impl IntoIterator<Item = &'p Point>) -> bool {
    pts.into_iter().all(|p| p.in_range())
}

/// Every vertex of a ring set.
pub fn ring_points(rings: &[Ring]) -> impl Iterator<Item = &Point> {
    rings.iter().flat_map(|r| r.iter())
}

/// Every vertex of a polygon set.
pub fn poly_points(polys: &[Polygon]) -> impl Iterator<Item = &Point> {
    polys.iter().flat_map(|p| p.rings().flat_map(|r| r.iter()))
}

/// Prints the decoded input to stderr when `POLYCLIP_FUZZ_DEBUG` is set (used to turn a
/// crash artifact into a readable test case: `POLYCLIP_FUZZ_DEBUG=1 cargo fuzz run
/// <target> <artifact>`).
pub fn dump(what: &str, v: &dyn core::fmt::Debug) {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    if *ON.get_or_init(|| std::env::var_os("POLYCLIP_FUZZ_DEBUG").is_some()) {
        eprintln!("{what} = {v:?};");
    }
}

/// Euclidean perimeter of all rings of a polygon set (f64).
pub fn perimeter(polys: &[Polygon]) -> f64 {
    polys
        .iter()
        .flat_map(|p| p.rings())
        .map(|r| {
            r.edges()
                .map(|(a, b)| ((b.x - a.x) as f64).hypot((b.y - a.y) as f64))
                .sum::<f64>()
        })
        .sum()
}
