//! Differential tests: polyclip vs Clipper2 on random and CAD-like inputs.
//!
//! Results are compared by total area and by the area of their symmetric difference
//! (computed exactly by polyclip). Both libraries round new vertices to the integer
//! grid, polyclip by snap rounding (every point moves by < 1 unit), Clipper2 by rounding
//! intersection points; the symmetric difference of two correct results is therefore a
//! band along the boundary whose area is bounded by a small multiple of the input
//! perimeter. Rectilinear inputs have integral intersections and must match exactly.
//!
//! Offsets additionally differ in how arcs are flattened (vertex counts and phases), so
//! their tolerance is `(arc tolerance + rounding) × perimeter`.
//!
//! Run with `cargo test --release -- --nocapture` to see per-configuration statistics.

use polyclip::*;
use polyclip_oracle::clipper;
use polyclip_oracle::compare::{Diff, as_rings, diff, perimeter, rings_of};
use polyclip_oracle::r#gen::*;

const RULES: [FillRule; 4] = [
    FillRule::EvenOdd,
    FillRule::NonZero,
    FillRule::Positive,
    FillRule::Negative,
];
const OPS: [Op; 4] = [Op::Union, Op::Intersection, Op::Difference, Op::Xor];

/// Per-configuration statistics; failures are collected and reported together.
#[derive(Default)]
struct Stats {
    cases: usize,
    max_ratio: f64,
    max_xor: f64,
    skipped: usize,
    failures: Vec<String>,
}

impl Stats {
    fn record(&mut self, d: Diff, what: impl FnOnce() -> String) {
        self.cases += 1;
        self.max_ratio = self.max_ratio.max(d.ratio());
        self.max_xor = self.max_xor.max(d.xor);
        if !d.ok() && self.failures.len() < 10 {
            self.failures.push(format!("{} -> {d:?}", what()));
        }
    }
    fn report(&self, name: &str) {
        eprintln!(
            "{name:<44} cases {:>4} (skipped {:>2})  max xor {:>16.1}  max xor/tol {:.4}  failures {}",
            self.cases,
            self.skipped,
            self.max_xor,
            self.max_ratio,
            self.failures.len()
        );
    }
}

fn finish(all: Vec<(String, Stats)>) {
    let mut failed = Vec::new();
    for (name, s) in &all {
        s.report(name);
        for f in &s.failures {
            failed.push(format!("[{name}] {f}"));
        }
    }
    assert!(
        failed.is_empty(),
        "{} mismatch(es) beyond tolerance:\n{}",
        failed.len(),
        failed.join("\n")
    );
}

fn ring_perimeter(rings: &[Ring]) -> f64 {
    perimeter(rings.iter().map(|r| r.as_slice()))
}

fn points(rings: &[Ring]) -> Vec<Vec<Point>> {
    rings.iter().map(|r| r.0.clone()).collect()
}

fn compare_boolean(
    stats: &mut Stats,
    a: &[Ring],
    b: &[Ring],
    rule: FillRule,
    op: Op,
    exact: bool,
    label: &dyn Fn() -> String,
) {
    let p = boolean(op, a, b, rule).expect("polyclip boolean");
    assert_eq!(check_canonical(&p, true), Ok(()), "{}", label());
    let c = clipper::boolean(op, &points(a), &points(b), rule);
    let tol = if exact {
        0.0
    } else {
        2.0 * (ring_perimeter(a) + ring_perimeter(b)) + 4.0
    };
    stats.record(diff(&p, &c, tol), || {
        format!("{} {op:?} {rule:?}\n  a = {a:?}\n  b = {b:?}", label())
    });
}

/// Random self-intersecting rings at several coordinate ranges, every op × fill rule.
#[test]
fn boolean_random_rings() {
    let mut all = Vec::new();
    for (range, cases) in [(8i64, 60), (1_000, 60), (1_000_000, 40), (MAX_COORD, 30)] {
        for rule in RULES {
            for op in OPS {
                let mut st = Stats::default();
                let mut s = Lcg(range as u64 ^ (rule as u64) << 8 ^ (op as u64) << 12);
                for case in 0..cases {
                    let seed = s.0;
                    let na = s.range(1, 3) as usize;
                    let nb = s.range(1, 3) as usize;
                    let a: Vec<Ring> = (0..na)
                        .map(|_| {
                            let n = s.range(3, 12) as usize;
                            random_ring(&mut s, n, range)
                        })
                        .collect();
                    let b: Vec<Ring> = (0..nb)
                        .map(|_| {
                            let n = s.range(3, 12) as usize;
                            random_ring(&mut s, n, range)
                        })
                        .collect();
                    compare_boolean(&mut st, &a, &b, rule, op, false, &|| {
                        format!("range {range} case {case} seed {seed:#x}")
                    });
                }
                all.push((format!("random ±{range} {op:?} {rule:?}"), st));
            }
        }
    }
    finish(all);
}

/// Rectangles only: every intersection is integral, so both libraries are exact.
#[test]
fn boolean_rectangles_exact() {
    let mut all = Vec::new();
    for rule in RULES {
        for op in OPS {
            let mut st = Stats::default();
            let mut s = Lcg(0x5eed ^ (rule as u64) << 8 ^ (op as u64) << 12);
            for case in 0..80 {
                let range = [10, 1000, 1 << 30][case % 3];
                let a: Vec<Ring> = (0..s.range(1, 8)).map(|_| rect(&mut s, range)).collect();
                let b: Vec<Ring> = (0..s.range(1, 8)).map(|_| rect(&mut s, range)).collect();
                compare_boolean(&mut st, &a, &b, rule, op, true, &|| format!("case {case}"));
            }
            all.push((format!("rects (exact) {op:?} {rule:?}"), st));
        }
    }
    finish(all);
}

/// CAD-like inputs (nm coordinates): rectangles, circles and tracks, many overlaps.
#[test]
fn boolean_cad_like() {
    let mut all = Vec::new();
    for rule in [FillRule::NonZero, FillRule::EvenOdd] {
        for op in OPS {
            let mut st = Stats::default();
            let mut s = Lcg(0xcad ^ (rule as u64) << 8 ^ (op as u64) << 12);
            for case in 0..40 {
                let range = 10_000_000; // ±10 mm
                let a: Vec<Ring> = (0..s.range(5, 30))
                    .map(|_| cad_ring(&mut s, range))
                    .collect();
                let b: Vec<Ring> = (0..s.range(5, 30))
                    .map(|_| cad_ring(&mut s, range))
                    .collect();
                compare_boolean(&mut st, &a, &b, rule, op, false, &|| format!("case {case}"));
            }
            all.push((format!("cad-like {op:?} {rule:?}"), st));
        }
    }
    finish(all);
}

fn compare_offset(
    stats: &mut Stats,
    input: &PolygonSet,
    delta: i64,
    join: Join,
    arc_tol: i64,
    label: &dyn Fn() -> String,
) {
    let p =
        offset(input, delta, join, ArcTol::new(arc_tol, Side::Inside)).expect("polyclip offset");
    assert_eq!(check_canonical(&p, true), Ok(()), "{}", label());
    let c = clipper::offset(&rings_of(input), delta as f64, join, arc_tol as f64);
    let per = perimeter(rings_of(&p).iter().map(|r| r.as_slice()))
        + perimeter(c.iter().map(|r| r.as_slice()));
    let tol = (arc_tol as f64 + 2.0) * per + 4.0;
    stats.record(diff(&p, &c, tol), || {
        format!(
            "{} delta {delta} {join:?} tol {arc_tol}\n  input = {input:?}",
            label()
        )
    });
}

/// Offsets of normalized random polygons (both signs) with round, square and bevel joins.
/// Miter joins are compared on CAD-like input only (below): the libraries clip
/// over-long miters differently (polyclip at `limit × delta`, Clipper2 squares at
/// `delta`), which is by design.
#[test]
fn offset_random_polygons() {
    let mut all = Vec::new();
    for join in [Join::Round, Join::Square, Join::Bevel] {
        for range in [1_000i64, 1_000_000] {
            let mut st = Stats::default();
            let mut s = Lcg(0x0ff5 ^ range as u64);
            for case in 0..60 {
                let rings: Vec<Ring> = (0..s.range(1, 3))
                    .map(|_| {
                        let n = s.range(3, 10) as usize;
                        random_ring(&mut s, n, range)
                    })
                    .collect();
                let input = union_all(&rings, FillRule::NonZero).unwrap();
                let delta = s.range(-range / 5, range / 5);
                let arc_tol = (range / 1000).max(1) * s.range(1, 5);
                compare_offset(&mut st, &input, delta, join, arc_tol, &|| {
                    format!("case {case}")
                });
            }
            all.push((format!("offset random ±{range} {join:?}"), st));
        }
    }
    finish(all);
}

/// Offsets of unions of CAD-like shapes, every join (miter limit 2: all convex corners of
/// these shapes are at most 90°, so both libraries miter every corner).
#[test]
fn offset_cad_like() {
    let mut all = Vec::new();
    for join in [
        Join::Round,
        Join::Square,
        Join::Bevel,
        Join::Miter { limit: 2.0 },
    ] {
        let mut st = Stats::default();
        let mut s = Lcg(0xcad0ff);
        for case in 0..40 {
            let range = 10_000_000;
            let rings: Vec<Ring> = (0..s.range(3, 20))
                .map(|_| cad_ring(&mut s, range))
                .collect();
            let input = union_all(&rings, FillRule::NonZero).unwrap();
            let delta = s.range(-500_000, 500_000);
            let arc_tol = s.range(100, 5_000);
            compare_offset(&mut st, &input, delta, join, arc_tol, &|| {
                format!("case {case}")
            });
        }
        all.push((format!("offset cad-like {join:?}"), st));
    }
    finish(all);
}

/// Every turn is gentle enough that a miter with limit 2 is never clipped: the miter
/// length is `delta / cos(turn / 2)`, at most `2 delta` for turns up to 120 degrees
/// (`cos(turn) >= -1/2`). Only then do the two libraries' miters coincide.
fn turns_within_miter_limit_2(pts: &[Point], closed: bool) -> bool {
    let mut v: Vec<Point> = pts.to_vec();
    v.dedup();
    if closed {
        while v.len() > 1 && v.first() == v.last() {
            v.pop();
        }
    }
    let n = v.len();
    if n < 3 {
        return !closed;
    }
    let dir = |a: Point, b: Point| {
        let (dx, dy) = ((b.x - a.x) as f64, (b.y - a.y) as f64);
        let l = dx.hypot(dy);
        (dx / l, dy / l)
    };
    let mut range = if closed { 0..n } else { 1..n - 1 };
    range.all(|i| {
        let (a, b, c) = (v[(i + n - 1) % n], v[i], v[(i + 1) % n]);
        let (d1, d2) = (dir(a, b), dir(b, c));
        d1.0 * d2.0 + d1.1 * d2.1 >= -0.49
    })
}

/// The closed loop of a path (for `EndCap::Joined`) is a simple ring.
fn simple_loop(pts: &[Point]) -> bool {
    let mut v: Vec<Point> = pts.to_vec();
    v.dedup();
    while v.len() > 1 && v.first() == v.last() {
        v.pop();
    }
    validate(&Polygon::new(Ring(v), vec![])).is_ok()
}

/// Strokes of open paths (random and rectilinear), every end cap and join.
///
/// Skipped (counted in the report), because the libraries legitimately differ there:
/// * miter joins at turns sharper than 120 degrees (over-long miters are clipped at
///   `limit × delta` by polyclip and squared at `delta` by Clipper2);
/// * `Joined` paths whose closed loop is not simple: Clipper2 strokes a closed path by
///   offsetting it outwards and inwards and filling with the positive rule, which drops
///   parts of the band around self-crossing loops, whereas polyclip returns the full
///   stroke (equal to the stroke of the explicitly closed open path).
#[test]
fn offset_paths_strokes() {
    let mut all = Vec::new();
    for cap in [EndCap::Round, EndCap::Square, EndCap::Butt, EndCap::Joined] {
        for join in [
            Join::Round,
            Join::Square,
            Join::Bevel,
            Join::Miter { limit: 2.0 },
        ] {
            let mut st = Stats::default();
            let mut s = Lcg(0x57 ^ (cap as u64) << 4);
            for case in 0..80 {
                let range = 1_000_000;
                let rectilinear = case % 2 == 0;
                let paths: Vec<Path> = (0..s.range(1, 3))
                    .map(|_| {
                        let n = s.range(2, 8) as usize;
                        if rectilinear {
                            rectilinear_path(&mut s, n, range)
                        } else {
                            random_path(&mut s, n, range)
                        }
                    })
                    .collect();
                let delta = s.range(1, range / 20);
                let arc_tol = s.range(10, 1000);
                let closed = cap == EndCap::Joined;
                if (closed && !paths.iter().all(|p| simple_loop(p)))
                    || (matches!(join, Join::Miter { .. })
                        && !paths.iter().all(|p| turns_within_miter_limit_2(p, closed)))
                {
                    st.skipped += 1;
                    continue;
                }
                let p = offset_paths(&paths, delta, join, cap, ArcTol::new(arc_tol, Side::Inside))
                    .expect("polyclip offset_paths");
                assert_eq!(check_canonical(&p, true), Ok(()));
                let pts: Vec<Vec<Point>> = paths.iter().map(|p| p.0.clone()).collect();
                let c = clipper::offset_paths(&pts, delta as f64, join, cap, arc_tol as f64);
                let per = perimeter(rings_of(&p).iter().map(|r| r.as_slice()))
                    + perimeter(c.iter().map(|r| r.as_slice()));
                let tol = (arc_tol as f64 + 2.0) * per + 4.0;
                st.record(diff(&p, &c, tol), || {
                    format!("case {case} delta {delta} {join:?} {cap:?} tol {arc_tol}\n  paths = {paths:?}")
                });
            }
            all.push((format!("stroke {cap:?} {join:?}"), st));
        }
    }
    finish(all);
}

/// Sanity check of the harness itself: Clipper2's orientation and fill-rule conventions
/// match polyclip's (counter-clockwise = positive winding, Y up).
#[test]
fn harness_conventions() {
    let ccw = Ring::from([(0, 0), (10, 0), (10, 10), (0, 10)]);
    let mut cw = ccw.clone();
    cw.reverse_orientation();
    for (r, pos) in [(&ccw, true), (&cw, false)] {
        let c = clipper::boolean(
            Op::Union,
            std::slice::from_ref(&r.0),
            &[],
            FillRule::Positive,
        );
        assert_eq!(c.is_empty(), !pos, "Positive rule on {r:?}");
        let p = union_all(r, FillRule::Positive).unwrap();
        assert_eq!(p.is_empty(), !pos);
    }
    let c = clipper::boolean(
        Op::Union,
        std::slice::from_ref(&ccw.0),
        &[],
        FillRule::NonZero,
    );
    assert!(diff(&union_all(&ccw, FillRule::NonZero).unwrap(), &c, 0.0).ok());
    assert_eq!(as_rings(&c).len(), 1);
    // Stroke half-width convention: a butt-capped segment of length 1000, delta 10, is a
    // 1000 x 20 rectangle in both.
    let seg = vec![vec![Point::new(0, 0), Point::new(1000, 0)]];
    let c = clipper::offset_paths(&seg, 10.0, Join::Round, EndCap::Butt, 1.0);
    assert_eq!(polyclip_oracle::compare::region_area(&c), 20_000.0);
}
