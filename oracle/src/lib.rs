//! Clipper2 as a differential-testing oracle for polyclip.
//!
//! * [`clipper`]: a minimal safe wrapper over the Clipper2 C API (`clipper2c-sys`): boolean
//!   operations, closed-polygon offsets and open-path offsets on `i64` coordinates.
//! * [`gen`]: deterministic input generators (random self-intersecting rings at several
//!   coordinate ranges, rectangles, circles, PCB-like tracks).
//! * [`compare`]: comparison by total area and by the area of the symmetric difference,
//!   computed exactly with polyclip itself, against a tolerance proportional to the
//!   perimeter.
//!
//! Nothing here is used by the polyclip library.

pub mod clipper {
    //! Safe wrapper over the parts of the Clipper2 C API the oracle needs.

    use clipper2c_sys::*;
    use polyclip::{EndCap, FillRule, Join, Op, Point};
    use std::os::raw::c_int;
    use std::time::{Duration, Instant};

    /// An owned `Paths64`.
    pub struct Paths(*mut ClipperPaths64);

    impl Paths {
        /// Empty path list.
        pub fn new() -> Self {
            // SAFETY: allocation of the documented size, then placement construction.
            unsafe {
                let mem = clipper_allocate(clipper_paths64_size());
                Paths(clipper_paths64(mem))
            }
        }

        /// Copies vertex lists into a new path list.
        pub fn from_points<'a>(paths: impl IntoIterator<Item = &'a [Point]>) -> Self {
            let out = Paths::new();
            for p in paths {
                let mut pts: Vec<ClipperPoint64> = p
                    .iter()
                    .map(|q| ClipperPoint64 { x: q.x, y: q.y })
                    .collect();
                // SAFETY: `pts` outlives the call; the C side copies the points into a new
                // Path64, which `add_path` copies again before we delete it.
                unsafe {
                    let mem = clipper_allocate(clipper_path64_size());
                    let path = clipper_path64_of_points(mem, pts.as_mut_ptr(), pts.len());
                    clipper_paths64_add_path(out.0, path);
                    clipper_delete_path64(path);
                }
            }
            out
        }

        /// Copies the paths back out.
        pub fn to_points(&self) -> Vec<Vec<Point>> {
            // SAFETY: indices stay within the lengths reported by the C side.
            unsafe {
                let n = clipper_paths64_length(self.0);
                (0..n)
                    .map(|i| {
                        let m = clipper_paths64_path_length(self.0, i as c_int);
                        (0..m)
                            .map(|j| {
                                let p = clipper_paths64_get_point(self.0, i as c_int, j as c_int);
                                Point::new(p.x, p.y)
                            })
                            .collect()
                    })
                    .collect()
            }
        }
    }

    impl Default for Paths {
        fn default() -> Self {
            Self::new()
        }
    }

    impl Drop for Paths {
        fn drop(&mut self) {
            // SAFETY: created by `clipper_paths64*` on `clipper_allocate` memory.
            unsafe { clipper_delete_paths64(self.0) }
        }
    }

    fn fill(rule: FillRule) -> ClipperFillRule {
        match rule {
            FillRule::EvenOdd => ClipperFillRule_EVEN_ODD,
            FillRule::NonZero => ClipperFillRule_NON_ZERO,
            FillRule::Positive => ClipperFillRule_POSITIVE,
            FillRule::Negative => ClipperFillRule_NEGATIVE,
        }
    }

    fn clip_type(op: Op) -> ClipperClipType {
        match op {
            Op::Union => ClipperClipType_UNION,
            Op::Intersection => ClipperClipType_INTERSECTION,
            Op::Difference => ClipperClipType_DIFFERENCE,
            Op::Xor => ClipperClipType_XOR,
        }
    }

    fn join_type(join: Join) -> (ClipperJoinType, f64) {
        match join {
            Join::Round => (ClipperJoinType_ROUND_JOIN, 2.0),
            Join::Bevel => (ClipperJoinType_BEVEL_JOIN, 2.0),
            Join::Square => (ClipperJoinType_SQUARE_JOIN, 2.0),
            Join::Miter { limit } => (ClipperJoinType_MITER_JOIN, limit),
        }
    }

    fn end_type(cap: EndCap) -> ClipperEndType {
        match cap {
            EndCap::Round => ClipperEndType_ROUND_END,
            EndCap::Square => ClipperEndType_SQUARE_END,
            EndCap::Butt => ClipperEndType_BUTT_END,
            EndCap::Joined => ClipperEndType_JOINED_END,
        }
    }

    /// `subject op clip` under `rule` (Clipper2 has one rule for both operands). Returns
    /// the closed output paths and the time spent in Clipper2's `Execute` alone.
    pub fn boolean_timed(
        op: Op,
        subject: &Paths,
        clip: &Paths,
        rule: FillRule,
    ) -> (Paths, Duration) {
        let closed = Paths::new();
        let open = Paths::new();
        // SAFETY: the engine is constructed, used and deleted here; the path lists are
        // valid for the duration (the engine copies them on `add_*`).
        let dt = unsafe {
            let mem = clipper_allocate(clipper_clipper64_size());
            let c = clipper_clipper64(mem);
            clipper_clipper64_add_subject(c, subject.0);
            clipper_clipper64_add_clip(c, clip.0);
            let t = Instant::now();
            let ok = clipper_clipper64_execute(c, clip_type(op), fill(rule), closed.0, open.0);
            let dt = t.elapsed();
            clipper_delete_clipper64(c);
            assert_eq!(ok, 1, "Clipper2 execute failed");
            dt
        };
        (closed, dt)
    }

    /// `subject op clip` under `rule`, as vertex lists.
    pub fn boolean(
        op: Op,
        subject: &[Vec<Point>],
        clip: &[Vec<Point>],
        rule: FillRule,
    ) -> Vec<Vec<Point>> {
        let s = Paths::from_points(subject.iter().map(|r| r.as_slice()));
        let c = Paths::from_points(clip.iter().map(|r| r.as_slice()));
        boolean_timed(op, &s, &c, rule).0.to_points()
    }

    fn offset_impl(
        paths: &Paths,
        delta: f64,
        join: Join,
        end: ClipperEndType,
        arc_tol: f64,
    ) -> (Paths, Duration) {
        let (jt, limit) = join_type(join);
        // SAFETY: as above; `execute` placement-constructs the result in fresh memory.
        unsafe {
            let mem = clipper_allocate(clipper_clipperoffset_size());
            let c = clipper_clipperoffset(mem, limit, arc_tol, 0, 0);
            clipper_clipperoffset_add_paths64(c, paths.0, jt, end);
            let out_mem = clipper_allocate(clipper_paths64_size());
            let t = Instant::now();
            let out = clipper_clipperoffset_execute(out_mem, c, delta);
            let dt = t.elapsed();
            let err = clipper_clipperoffset_error_code(c);
            clipper_delete_clipperoffset(c);
            assert_eq!(err, 0, "Clipper2 offset error {err}");
            (Paths(out), dt)
        }
    }

    /// Offsets closed polygons (outer rings counter-clockwise, holes clockwise) by
    /// `delta`, with the given join, miter limit (from `Join::Miter`) and arc tolerance (in
    /// coordinate units). Returns the result and the time spent in `Execute`.
    pub fn offset_timed(rings: &Paths, delta: f64, join: Join, arc_tol: f64) -> (Paths, Duration) {
        offset_impl(rings, delta, join, ClipperEndType_POLYGON_END, arc_tol)
    }

    /// Closed-polygon offset as vertex lists. See [`offset_timed`].
    pub fn offset(rings: &[Vec<Point>], delta: f64, join: Join, arc_tol: f64) -> Vec<Vec<Point>> {
        let p = Paths::from_points(rings.iter().map(|r| r.as_slice()));
        offset_timed(&p, delta, join, arc_tol).0.to_points()
    }

    /// Open-path offset (stroke of half-width `delta`) as vertex lists.
    pub fn offset_paths(
        paths: &[Vec<Point>],
        delta: f64,
        join: Join,
        cap: EndCap,
        arc_tol: f64,
    ) -> Vec<Vec<Point>> {
        let p = Paths::from_points(paths.iter().map(|r| r.as_slice()));
        offset_impl(&p, delta, join, end_type(cap), arc_tol)
            .0
            .to_points()
    }
}

pub mod r#gen {
    //! Deterministic input generators.

    use polyclip::{Path, Point, Ring};
    use std::f64::consts::PI;

    /// Small, fast, deterministic PRNG (64-bit LCG, high bits).
    #[derive(Clone, Debug)]
    pub struct Lcg(pub u64);

    impl Lcg {
        /// Next 31 random bits.
        pub fn bits(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0 >> 33
        }
        /// Uniform in `lo..=hi`.
        pub fn range(&mut self, lo: i64, hi: i64) -> i64 {
            let span = (hi - lo) as u128 + 1;
            let r = ((self.bits() as u128) << 31 | self.bits() as u128) % span;
            lo + r as i64
        }
        /// Uniform in `[0, 1)`.
        pub fn unit(&mut self) -> f64 {
            self.bits() as f64 / (1u64 << 31) as f64
        }
        /// A random boolean.
        pub fn coin(&mut self) -> bool {
            self.bits() & 1 == 1
        }
    }

    /// Ring with `n` random vertices in `[-range, range]^2` (usually self-intersecting).
    pub fn random_ring(s: &mut Lcg, n: usize, range: i64) -> Ring {
        (0..n)
            .map(|_| Point::new(s.range(-range, range), s.range(-range, range)))
            .collect()
    }

    /// Axis-aligned rectangle inside `[-range, range]^2`, random orientation.
    pub fn rect(s: &mut Lcg, range: i64) -> Ring {
        let (x0, y0) = (s.range(-range, range), s.range(-range, range));
        let (w, h) = (s.range(1, range), s.range(1, range));
        let mut r = Ring::from([(x0, y0), (x0 + w, y0), (x0 + w, y0 + h), (x0, y0 + h)]);
        if s.coin() {
            r.reverse_orientation();
        }
        r
    }

    /// Circle approximated by `n` rounded vertices (counter-clockwise).
    pub fn circle(cx: i64, cy: i64, r: f64, n: usize) -> Ring {
        (0..n)
            .map(|k| {
                let a = 2.0 * PI * k as f64 / n as f64;
                Point::new(
                    cx + (r * a.cos()).round() as i64,
                    cy + (r * a.sin()).round() as i64,
                )
            })
            .collect()
    }

    /// Random circle inside `[-range, range]^2`.
    pub fn random_circle(s: &mut Lcg, range: i64) -> Ring {
        let r = (s.range(range / 50 + 1, range / 4 + 2)) as f64;
        let n = [16, 32, 64][s.range(0, 2) as usize];
        circle(s.range(-range, range), s.range(-range, range), r, n)
    }

    /// PCB-like track: a stadium (segment with round caps) of half-width `w`.
    pub fn track(a: Point, b: Point, w: f64, cap_pts: usize) -> Ring {
        let (dx, dy) = ((b.x - a.x) as f64, (b.y - a.y) as f64);
        let a0 = dy.atan2(dx) - PI / 2.0;
        let mut v = Vec::with_capacity(2 * cap_pts + 2);
        for (c, start) in [(b, a0), (a, a0 + PI)] {
            for k in 0..=cap_pts {
                let t = start + PI * k as f64 / cap_pts as f64;
                v.push(Point::new(
                    c.x + (w * t.cos()).round() as i64,
                    c.y + (w * t.sin()).round() as i64,
                ));
            }
        }
        v.dedup();
        Ring(v)
    }

    /// Random track inside `[-range, range]^2`, mostly at 0/45/90 degrees like real
    /// routing.
    pub fn random_track(s: &mut Lcg, range: i64) -> Ring {
        let a = Point::new(s.range(-range, range), s.range(-range, range));
        let len = s.range(range / 20 + 1, range / 2 + 2);
        let (dx, dy) = match s.range(0, 4) {
            0 => (len, 0),
            1 => (0, len),
            2 => (len, len),
            3 => (len, -len),
            _ => (s.range(-len, len), s.range(-len, len)),
        };
        let b = Point::new(a.x + dx, a.y + dy);
        let w = s.range(range / 200 + 1, range / 40 + 2) as f64;
        track(a, b, w, 8)
    }

    /// A CAD-like ring: rectangle, circle or track.
    pub fn cad_ring(s: &mut Lcg, range: i64) -> Ring {
        match s.range(0, 2) {
            0 => rect(s, range),
            1 => random_circle(s, range),
            _ => random_track(s, range),
        }
    }

    /// Random polyline of `n` vertices.
    pub fn random_path(s: &mut Lcg, n: usize, range: i64) -> Path {
        Path(random_ring(s, n, range).0)
    }

    /// Rectilinear polyline (all 90-degree turns).
    pub fn rectilinear_path(s: &mut Lcg, n: usize, range: i64) -> Path {
        let mut p = vec![Point::new(s.range(-range, range), s.range(-range, range))];
        for i in 1..n {
            let l = p[i - 1];
            let d = s.range(range / 10 + 1, range / 2 + 2) * if s.coin() { 1 } else { -1 };
            p.push(if i % 2 == 0 {
                Point::new(l.x + d, l.y)
            } else {
                Point::new(l.x, l.y + d)
            });
        }
        Path(p)
    }
}

pub mod compare {
    //! Area-based comparison of two results.

    use polyclip::{Boolean, FillRule, Op, Point, PolygonSet, Ring, area2};

    /// Rings of a polygon set as plain vertex lists.
    pub fn rings_of(set: &PolygonSet) -> Vec<Vec<Point>> {
        set.iter()
            .flat_map(|p| p.rings().map(|r| r.0.clone()))
            .collect()
    }

    /// Converts vertex lists to rings.
    pub fn as_rings(v: &[Vec<Point>]) -> Vec<Ring> {
        v.iter().map(|r| Ring(r.clone())).collect()
    }

    /// Euclidean length of the closed rings.
    pub fn perimeter<'a>(rings: impl IntoIterator<Item = &'a [Point]>) -> f64 {
        rings
            .into_iter()
            .map(|r| {
                let n = r.len();
                (0..n)
                    .map(|i| {
                        let (a, b) = (r[i], r[(i + 1) % n]);
                        ((b.x - a.x) as f64).hypot((b.y - a.y) as f64)
                    })
                    .sum::<f64>()
            })
            .sum()
    }

    /// Area of a region given as rings under the non-zero rule (any orientation), exact.
    pub fn region_area(rings: &[Vec<Point>]) -> f64 {
        let set = polyclip::union_all(&as_rings(rings), FillRule::NonZero).expect("in range");
        area2(&set) as f64 / 2.0
    }

    /// Outcome of comparing a polyclip result with a Clipper2 result.
    #[derive(Clone, Copy, Debug)]
    pub struct Diff {
        /// polyclip area.
        pub area_p: f64,
        /// Clipper2 area (normalized with the non-zero rule).
        pub area_c: f64,
        /// Area of the symmetric difference.
        pub xor: f64,
        /// Allowed deviation.
        pub tol: f64,
    }

    impl Diff {
        /// Both the area difference and the symmetric-difference area are within `tol`.
        pub fn ok(&self) -> bool {
            (self.area_p - self.area_c).abs() <= self.tol && self.xor <= self.tol
        }
        /// `xor / tol` (how much of the tolerance is used).
        pub fn ratio(&self) -> f64 {
            if self.tol > 0.0 {
                self.xor / self.tol
            } else {
                self.xor
            }
        }
    }

    /// Compares `p` (polyclip) with `c` (Clipper2 paths, outer rings positive, holes
    /// negative) by area and symmetric-difference area, with tolerance `tol` (area units).
    pub fn diff(p: &PolygonSet, c: &[Vec<Point>], tol: f64) -> Diff {
        let cr = as_rings(c);
        let x = Boolean::new()
            .subject(p, FillRule::NonZero)
            .clip(&cr, FillRule::NonZero)
            .op(Op::Xor)
            .execute()
            .expect("in range");
        Diff {
            area_p: area2(p) as f64 / 2.0,
            area_c: region_area(c),
            xor: area2(&x) as f64 / 2.0,
            tol,
        }
    }
}
