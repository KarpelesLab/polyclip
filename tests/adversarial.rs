//! Adversarial regression tests.
//!
//! Tests marked `#[ignore = "bug: ..."]` reproduce confirmed bugs (wrong results or
//! performance / memory cliffs) and fail today; run them with
//! `cargo test --test adversarial -- --ignored`. Remove the `#[ignore]` once fixed.
//! The other tests are passing regression coverage for adversarial inputs that were
//! investigated and found to be handled correctly.

use polyclip::*;
use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;
use std::time::{Duration, Instant};

fn p(x: i64, y: i64) -> Point {
    Point::new(x, y)
}

fn rect(x0: i64, y0: i64, w: i64, h: i64) -> Ring {
    Ring::from([(x0, y0), (x0 + w, y0), (x0 + w, y0 + h), (x0, y0 + h)])
}

/// Deterministic LCG for the self-contained randomized checks.
struct Rng(u64);

impl Rng {
    fn r(&mut self, m: i64) -> i64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 33) as i64).rem_euclid(m)
    }
}

// ---------------------------------------------------------------------------------------
// Allocation tracking: the largest single allocation made by the current thread, so a
// memory blow-up can be asserted deterministically (tests run on separate threads).

struct MaxAlloc;

thread_local! {
    static MAX_ALLOC: Cell<usize> = const { Cell::new(0) };
}

fn note_alloc(size: usize) {
    let _ = MAX_ALLOC.try_with(|m| {
        if size > m.get() {
            m.set(size)
        }
    });
}

// SAFETY: forwards every call to the system allocator unchanged; only records sizes.
unsafe impl GlobalAlloc for MaxAlloc {
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        note_alloc(l.size());
        // SAFETY: same contract as the caller's.
        unsafe { System.alloc(l) }
    }
    unsafe fn alloc_zeroed(&self, l: Layout) -> *mut u8 {
        note_alloc(l.size());
        // SAFETY: same contract as the caller's.
        unsafe { System.alloc_zeroed(l) }
    }
    unsafe fn dealloc(&self, ptr: *mut u8, l: Layout) {
        // SAFETY: same contract as the caller's.
        unsafe { System.dealloc(ptr, l) }
    }
    unsafe fn realloc(&self, ptr: *mut u8, l: Layout, new_size: usize) -> *mut u8 {
        note_alloc(new_size);
        // SAFETY: same contract as the caller's.
        unsafe { System.realloc(ptr, l, new_size) }
    }
}

#[global_allocator]
static ALLOC: MaxAlloc = MaxAlloc;

/// Largest single allocation (bytes) made by `f` on this thread.
fn max_alloc_during(f: impl FnOnce()) -> usize {
    MAX_ALLOC.with(|m| m.set(0));
    f();
    MAX_ALLOC.with(|m| m.get())
}

// =======================================================================================
// Confirmed bugs
// =======================================================================================

/// Expected: with `Side::Nearest` the polyline stays within `tolerance` of the true arc on
/// both sides. Actual: the arc's first and last segments (or the single chord when the arc
/// gets only one segment) cut into the disk by up to ~2x the tolerance (18.9 for
/// tolerance 10 here).
///
/// Root cause: src/arc.rs:471-480 (`Construction::Mid`) puts the interior vertices at
/// radius `rv = 2r / (1 + cos(h/2))`, but the arc's end points stay on the true circle.
/// `step_for` (src/arc.rs:356-361, `c = (r - t) / (r + t)`) sizes the step assuming every
/// vertex is at `rv`, so the end chords sag `r (1 - c) ≈ 2t` inside the circle.
#[test]
#[ignore = "bug: Side::Nearest arcs deviate inward by up to ~2x the tolerance"]
fn nearest_arc_chord_exceeds_tolerance() {
    // Pie slice, radius 1000, sweep 0.3896 rad < step 2*acos(990/1010) = 0.3987 rad: the
    // arc becomes the single chord (1000,0)-(925,380), sagitta ~18.9.
    let shape = Shape::new(
        vec![
            Curve::Line(p(1000, 0)),
            Curve::CenterArc {
                center: p(0, 0),
                end: p(925, 380),
                ccw: true,
            },
            Curve::Line(p(0, 0)),
        ],
        vec![],
    );
    let poly = shape.to_polygon(ArcTol::new(10, Side::Nearest)).unwrap();
    // (981, 194) lies on the true arc (radius 999.998).
    let boundary: Path = poly.outer.clone().into();
    let d = distance_sq(&p(981, 194), &boundary).unwrap().distance_f64();
    assert!(d <= 10.0 + 1.0, "deviation {d} exceeds tolerance 10"); // actual 18.92
}

/// Same defect through offsetting: round joins with `Side::Nearest` leave points of the
/// true offset arc up to ~1.5x the tolerance outside the result (15.15 for tolerance 10).
/// Root cause as in [`nearest_arc_chord_exceeds_tolerance`].
#[test]
#[ignore = "bug: Side::Nearest round joins deviate inward by more than the tolerance"]
fn nearest_round_join_exceeds_tolerance() {
    let sq = rect(0, 0, 10_000, 10_000);
    let g = offset(&sq, 1000, Join::Round, ArcTol::new(10, Side::Nearest)).unwrap();
    let boundary: Path = g[0].outer.clone().into();
    let mut worst = 0.0f64;
    for k in 0..=1000 {
        let a = std::f64::consts::FRAC_PI_2 * k as f64 / 1000.0;
        let q = p(
            10_000 + (1000.0 * a.cos()).round() as i64,
            10_000 + (1000.0 * a.sin()).round() as i64,
        );
        if locate(&g, q) == Location::Outside {
            worst = worst.max(distance_sq(&q, &boundary).unwrap().distance_f64());
        }
    }
    assert!(worst <= 11.0, "true-arc point {worst} outside the result"); // actual 15.15
}

/// Expected: a nearly flat three-point arc whose whole sagitta (1) is within the tolerance
/// (1) becomes a single chord. Actual: `Err(InvalidParameter("arc tolerance too small for
/// the radius"))` for every `Side` (radius ~1.25e17 for a 1 m chord with the mid point
/// 1 nm off the line; plausible with noisy imported arc data).
///
/// Root cause: `step_for` (src/arc.rs:356-367) computes `1 - t/r`, `r/(r+t)` or
/// `(r-t)/(r+t)`, which round to exactly 1.0 once `t/r` is below ~1.1e-16, so `acos`
/// gives a zero step and `segment_count` (src/arc.rs:396) sees an infinite count.
#[test]
#[ignore = "bug: nearly flat 3-point arc (huge radius) errors instead of becoming a chord"]
fn nearly_flat_arc_spurious_error() {
    let l = 500_000_000i64;
    let shape = Shape::new(
        vec![
            Curve::Arc {
                mid: p(0, 1),
                end: p(l, 0),
            },
            Curve::Line(p(l, -1000)),
            Curve::Line(p(-l, -1000)),
            Curve::Line(p(-l, 0)),
        ],
        vec![],
    );
    for side in [Side::Inside, Side::Outside, Side::Nearest] {
        assert!(shape.to_polygon(ArcTol::new(1, side)).is_ok(), "{side:?}");
    }
}

/// Expected: `minkowski_sum` of a self-intersecting ring equals the sum of the region it
/// denotes (non-zero rule; the general path normalizes with `union_all`). Actual: a
/// pentagram (convex pentagon vertices visited every other one) turns left at every
/// vertex, so `is_convex` (src/hull.rs:58-65, which only checks that every turn agrees in
/// sign with the signed area) accepts it and the convex shortcut (src/hull.rs:97-105)
/// returns the hull of the pairwise vertex sums: doubled area 49098 instead of 24256.
#[test]
#[ignore = "bug: minkowski_sum treats a self-intersecting star ring as convex"]
fn minkowski_sum_pentagram_is_not_convex() {
    let star = Polygon::from(Ring::from([
        (0, 100),
        (-59, -81),
        (95, 31),
        (-95, 31),
        (59, -81),
    ]));
    let square = Polygon::from(rect(-1, -1, 2, 2));
    let region = union_all(&star, FillRule::NonZero).unwrap();
    assert_eq!(region.len(), 1);
    let expected = minkowski_sum(&region[0], &square).unwrap();
    let got = minkowski_sum(&star, &square).unwrap();
    assert_eq!(area2(&got), area2(&expected));
    assert_eq!(got, expected);
}

/// Memory cliff in the noding index. Two tiny triangles at opposite ends of the
/// coordinate range on one row (bounding box 2^41 x 2) make `snap_round` allocate ~160 MB
/// in a single block; 2000 such triangles (6000 edges) peak at ~3.5 GB RSS, and a few
/// thousand more abort on allocation failure.
///
/// Root cause: `Grid::build_uniform` (src/node.rs:100-102) bounds the cell size from below
/// by `sqrt(w * h / (n + 16))`, which limits the *area* per cell but not the cell count
/// when the bounding box is very thin: with `h` ~ 1 the grid has
/// `nx ≈ sqrt(w * (n + 16) / h)` cells (here ~5e6 for 6 edges, ~8e7 for 6000), each with a
/// 32-byte leaf rectangle and CSR slots. The budget check only counts cells *listed by
/// segments*, not the cells allocated. The cell size should also be at least
/// `max(w, h) / (n + 16)` (or the grid rejected when `nx * ny` greatly exceeds `n`).
#[test]
#[ignore = "bug: snap_round allocates O(sqrt(width * n / height)) grid cells for thin inputs"]
fn thin_far_apart_input_allocates_huge_grid() {
    let m = MAX_COORD;
    let a = vec![
        Ring::from([(-m, 0), (-m + 10, 0), (-m + 10, 1)]),
        Ring::from([(m - 10, 0), (m, 0), (m, 1)]),
    ];
    let mut r = Vec::new();
    let biggest = max_alloc_during(|| r = union_all(&a, FillRule::NonZero).unwrap());
    assert_eq!(r.len(), 2);
    assert!(
        biggest < 16 << 20,
        "largest single allocation {} MB for a 6-edge input",
        biggest >> 20
    );
}

/// Quadratic noding for parallel segments at a "generic" angle. 4000 parallel tracks of
/// slope 1/3 (thin parallelograms, no crossings at all) take ~0.9 s in release (16000:
/// 14 s), while the same bus at 45 degrees takes 5 ms (16000: 21 ms). Plausible CAD input:
/// a rotated bus or hatch lines at an arbitrary angle (`clip_paths` of slope-1/2 hatch
/// lines over a zone shows the same O(n^2) growth).
///
/// Root cause: the uniform grid is rejected (long segments exceed the cell budget) and the
/// k-d partition cannot split long parallel segments that all straddle every split line
/// (src/node.rs:302-305, "no progress"), so they share one leaf; `for_each_pair`
/// (src/node.rs:375-412) then sweeps along the best of four fixed directions (x, y, both
/// diagonals), none of which separates segments that are not aligned with them, and tests
/// all O(n^2) bounding-box-overlapping pairs. Radial patterns hit the same limit: a fan of
/// 8000 / 32000 pie slices with a common apex takes 1.2 s / 19 s to union in release.
#[test]
#[ignore = "bug: noding is quadratic for long parallel segments not aligned to x/y/diagonals"]
fn skewed_parallel_bus_is_quadratic() {
    let n = 4000i64;
    let bus = |dy: i64| -> Vec<Ring> {
        (0..n)
            .map(|i| {
                let y = i * 100;
                Ring::from([
                    (0, y),
                    (30_000_000, y + dy),
                    (30_000_000, y + dy + 40),
                    (0, y + 40),
                ])
            })
            .collect()
    };
    // Reference: the 45 degree bus (fast).
    let t0 = Instant::now();
    assert_eq!(
        union_all(&bus(30_000_000), FillRule::NonZero)
            .unwrap()
            .len(),
        n as usize
    );
    let diag = t0.elapsed();
    let t0 = Instant::now();
    assert_eq!(
        union_all(&bus(10_000_000), FillRule::NonZero)
            .unwrap()
            .len(),
        n as usize
    );
    let skew = t0.elapsed();
    assert!(
        skew < diag * 20 + Duration::from_millis(200),
        "slope 1/3: {skew:?}, slope 1: {diag:?}"
    );
}

/// `trapezoids` re-sorts the whole sweep status and re-matches every open trapezoid with a
/// linear `open.iter().position(..)` in every slab (src/decompose.rs:107-135). So n stacked
/// horizontal stripes with distinct x coordinates (2n slabs, ~2n active edges, n open
/// trapezoids) cost O(n^3), although the output has only n trapezoids. Release: 800
/// stripes 33 ms, 1600 0.23 s, 3200 1.6 s (x7 per doubling); debug, 3000 stripes: 14.8 s.
#[test]
#[ignore = "bug: trapezoids is cubic in the number of stacked parallel stripes"]
fn trapezoids_parallel_stripes_cubic() {
    let n = 3000i64;
    let rings: Vec<Ring> = (0..n)
        .map(|i| {
            Ring::from([
                (i, 10 * i),
                (1_000_000 + i, 10 * i),
                (1_000_000 + i, 10 * i + 5),
                (i, 10 * i + 5),
            ])
        })
        .collect();
    let t0 = Instant::now();
    let t = trapezoids(&rings).unwrap();
    assert_eq!(t.len(), n as usize);
    assert!(t0.elapsed() < Duration::from_secs(2), "{:?}", t0.elapsed());
}

/// `intersects` and `distance_less_than` find candidate pairs with `any_pair`
/// (src/query.rs:406-460), a sweep-and-prune along x only. Long horizontal edges all stay
/// in the active lists, so two interleaved sets of n parallel horizontal stripes (a "two
/// layers of parallel tracks" DRC query) test all O(n^2) pairs. Release, each call: 2000
/// stripes per side 0.11 s, 8000 1.7 s, 32000 27 s; debug, 8000: 31.8 s.
#[test]
#[ignore = "bug: intersects/distance_less_than are quadratic on parallel horizontal stripes"]
fn intersects_parallel_stripes_quadratic() {
    let n = 8000i64;
    let stripe = |y: i64| Polygon::from(rect(0, y, 1_000_000, 5));
    let a: Vec<Polygon> = (0..n).map(|i| stripe(20 * i)).collect();
    let b: Vec<Polygon> = (0..n).map(|i| stripe(20 * i + 10)).collect();
    let t0 = Instant::now();
    assert!(!intersects(&a, &b));
    assert!(!distance_less_than(&a, &b, 5));
    assert!(t0.elapsed() < Duration::from_secs(2), "{:?}", t0.elapsed());
}

/// `distance` prunes candidate pairs only by x (src/distance.rs:242-262: `sb` sorted by
/// min x, the inner loop breaks on `y.min.x > x.max.x + r`, `sa` unsorted, no y pruning).
/// For two long paths separated vertically by more than their segment length, every
/// segment of `a` is paired with every segment of `b` to its left: O(n*m). Release: 5000
/// vertices 14 ms, 20000 0.21 s, 80000 3.3 s, while `distance_less_than` on the same input
/// takes 1.3 ms; debug, 80000: 47 s.
#[test]
#[ignore = "bug: distance is quadratic for vertically separated long paths"]
fn distance_parallel_paths_quadratic() {
    let n = 80_000i64;
    let a: Path = (0..n).map(|i| p(i * 100, (i % 2) * 50)).collect();
    let b: Path = (0..n).map(|i| p(i * 100, 10_000 + (i % 2) * 50)).collect();
    let t0 = Instant::now();
    let c = distance(&a, &b).unwrap();
    assert_eq!(c.sq.cmp_dist(9_950), core::cmp::Ordering::Greater);
    assert!(t0.elapsed() < Duration::from_secs(2), "{:?}", t0.elapsed());
}

// =======================================================================================
// Passing regression coverage
// =======================================================================================

fn translate(ps: &PolygonSet, dx: i64, dy: i64) -> PolygonSet {
    let t = |r: &Ring| -> Ring { r.iter().map(|q| p(q.x + dx, q.y + dy)).collect() };
    ps.iter()
        .map(|poly| Polygon {
            outer: t(&poly.outer),
            holes: poly.holes.iter().map(t).collect(),
        })
        .collect()
}

fn random_walk_rings(rng: &mut Rng, n: usize, pts: usize, range: i64, step: i64) -> Vec<Ring> {
    (0..n)
        .map(|_| {
            let (mut x, mut y) = (rng.r(range), rng.r(range));
            (0..pts)
                .map(|_| {
                    x = (x + rng.r(2 * step + 1) - step).clamp(0, range);
                    y = (y + rng.r(2 * step + 1) - step).clamp(0, range);
                    p(x, y)
                })
                .collect()
        })
        .collect()
}

/// Snap rounding on the integer grid is translation invariant, and the result must not
/// depend on the spatial index (perturbed here by adding a far-away copy, which changes
/// the grid / k-d structure), on ring order or on ring start vertices. Checked on heavily
/// self-intersecting random walks, including coordinates near `MAX_COORD`.
#[test]
fn union_is_translation_index_and_order_invariant() {
    let mut rng = Rng(1);
    for it in 0..60 {
        let range: i64 = [50, 300, 3000, (1 << 40) - 1000][it % 4];
        let step = [3, 20, 200, 1 << 38][rng.r(4) as usize].min(range);
        let (nr, np) = (1 + rng.r(12) as usize, 3 + rng.r(40) as usize);
        let a = random_walk_rings(&mut rng, nr, np, range, step);
        let rule = [FillRule::NonZero, FillRule::EvenOdd, FillRule::Positive][rng.r(3) as usize];
        let r = union_all(&a, rule).unwrap();
        assert_eq!(check_canonical(&r, true), Ok(()), "it {it}");
        let (dx, dy) = (rng.r(1000) - 500, rng.r(1000) - 500);
        let moved: Vec<Ring> = a
            .iter()
            .map(|r| r.iter().map(|q| p(q.x + dx, q.y + dy)).collect())
            .collect();
        assert_eq!(
            union_all(&moved, rule).unwrap(),
            translate(&r, dx, dy),
            "it {it}"
        );
        if range < 1_000_000 {
            let far = 10_000_000;
            let mut b = a.clone();
            b.extend(a.iter().map(|r| {
                r.iter()
                    .map(|q| p(q.x + far, q.y + far / 3))
                    .collect::<Ring>()
            }));
            let mut exp = r.clone();
            exp.extend(translate(&r, far, far / 3));
            assert_eq!(union_all(&b, rule).unwrap(), exp, "it {it}");
        }
        let mut c: Vec<Ring> = a.iter().rev().cloned().collect();
        for (k, rr) in c.iter_mut().enumerate() {
            let l = rr.len();
            rr.rotate_left(k % l);
        }
        assert_eq!(union_all(&c, rule).unwrap(), r, "it {it}");
    }
}

/// `clip_paths` with 45-degree and axis-parallel paths against rectangles, all on even
/// coordinates: every crossing is an integer point and every segment meets only pixels
/// whose centres lie on it, so nothing is rerouted and the classification is exact. Checks
/// inside/outside classification of every piece edge, total length conservation and that
/// every piece edge carries the tag of a source edge containing it.
#[test]
fn clip_paths_diagonal_exact() {
    let mut rng = Rng(7);
    for it in 0..400 {
        let clip: Vec<Ring> = (0..rng.r(5))
            .map(|_| {
                let mut r = rect(
                    2 * rng.r(20),
                    2 * rng.r(20),
                    2 + 2 * rng.r(10),
                    2 + 2 * rng.r(10),
                );
                if rng.r(2) == 0 {
                    r.reverse_orientation();
                }
                r
            })
            .collect();
        let rule = [
            FillRule::NonZero,
            FillRule::EvenOdd,
            FillRule::Positive,
            FillRule::Negative,
        ][rng.r(4) as usize];
        let paths: Vec<TaggedPath> = (0..1 + rng.r(4))
            .map(|_| {
                let mut v = vec![p(2 * rng.r(20) - 4, 2 * rng.r(20) - 4)];
                for _ in 0..1 + rng.r(6) {
                    let l = *v.last().unwrap();
                    let d = 2 * (rng.r(13) - 6);
                    v.push(match rng.r(4) {
                        0 => p(l.x + d, l.y),
                        1 => p(l.x, l.y + d),
                        2 => p(l.x + d, l.y + d),
                        _ => p(l.x + d, l.y - d),
                    });
                }
                let tags = (1..v.len()).map(|_| rng.r(3) as u64).collect();
                TaggedPath { points: v, tags }
            })
            .collect();
        let r = clip_paths(&paths, &clip, rule).unwrap();
        let dbl: PolygonSet = union_all(&clip, rule)
            .unwrap()
            .iter()
            .map(|poly| {
                let d = |r: &Ring| -> Ring { r.iter().map(|q| p(2 * q.x, 2 * q.y)).collect() };
                Polygon {
                    outer: d(&poly.outer),
                    holes: poly.holes.iter().map(d).collect(),
                }
            })
            .collect();
        let l1 = |a: Point, b: Point| (b.x - a.x).abs() + (b.y - a.y).abs();
        let total_in: i64 = paths
            .iter()
            .map(|pp| pp.points.windows(2).map(|w| l1(w[0], w[1])).sum::<i64>())
            .sum();
        let mut total_out = 0;
        for (pieces, inside) in [(&r.inside, true), (&r.outside, false)] {
            for pc in pieces {
                for (k, w) in pc.points.windows(2).enumerate() {
                    total_out += l1(w[0], w[1]);
                    let loc = locate(&dbl, p(w[0].x + w[1].x, w[0].y + w[1].y));
                    assert_eq!(loc != Location::Outside, inside, "it {it}: {w:?}");
                    let t = pc.tags[k];
                    assert!(
                        paths.iter().any(|pp| {
                            pp.points.windows(2).zip(&pp.tags).any(|(s, &st)| {
                                st == t
                                    && predicates::on_segment(s[0], s[1], w[0])
                                    && predicates::on_segment(s[0], s[1], w[1])
                            })
                        }),
                        "it {it}: tag {t} on {w:?}"
                    );
                }
            }
        }
        assert_eq!(total_in, total_out, "it {it}");
    }
}

/// The nesting tree from `execute_tree` is consistent with geometry: every ring is
/// strictly inside exactly as many other rings as it has ancestors, holes sit at odd depth,
/// and parent / child links agree (`check_canonical` cannot see tree errors, since an island
/// wrongly attached as a root yields the same polygon set). Also checks tag determinism
/// under input reordering.
#[test]
fn polytree_nesting_matches_geometry() {
    let mut rng = Rng(3);
    for it in 0..1500 {
        let range = [12, 40, 1000][it % 3];
        let mut mk = || -> Vec<TaggedRing> {
            (0..rng.r(8))
                .map(|_| {
                    let pts: Vec<Point> = if rng.r(2) == 0 {
                        let r = rect(
                            rng.r(range),
                            rng.r(range),
                            1 + rng.r(range / 2),
                            1 + rng.r(range / 2),
                        );
                        r.0
                    } else {
                        (0..3 + rng.r(6))
                            .map(|_| p(rng.r(range), rng.r(range)))
                            .collect()
                    };
                    let tags = pts.iter().map(|_| rng.r(4) as u64).collect();
                    TaggedRing { points: pts, tags }
                })
                .collect()
        };
        let (a, b) = (mk(), mk());
        let op = [Op::Union, Op::Intersection, Op::Difference, Op::Xor][rng.r(4) as usize];
        let rules = [
            FillRule::NonZero,
            FillRule::EvenOdd,
            FillRule::Positive,
            FillRule::Negative,
        ];
        let (ra, rb) = (rules[rng.r(4) as usize], rules[rng.r(4) as usize]);
        let keep = rng.r(4) == 0;
        let run = |a: &[TaggedRing], b: &[TaggedRing]| {
            Boolean::new()
                .subject(a, ra)
                .clip(b, rb)
                .op(op)
                .keep_collinear(keep)
                .execute_tree()
                .unwrap()
        };
        let t = run(&a, &b);
        assert_eq!(
            check_canonical(&t.to_polygon_set(), false),
            Ok(()),
            "it {it}"
        );
        let n = t.nodes.len();
        for i in 0..n {
            let ring = &t.nodes[i].ring;
            let m = p(ring[0].x + ring[1].x, ring[0].y + ring[1].y);
            let mut inside = 0;
            for (j, other) in t.nodes.iter().enumerate() {
                if j == i {
                    continue;
                }
                let d: Vec<Point> = other.ring.iter().map(|q| p(2 * q.x, 2 * q.y)).collect();
                match locate_in_ring(&d, m) {
                    Location::Inside => inside += 1,
                    Location::OnBoundary => panic!("it {it}: rings {i} and {j} share an edge"),
                    Location::Outside => {}
                }
            }
            let (mut depth, mut cur, mut child) = (0, t.nodes[i].parent, i);
            while let Some(par) = cur {
                assert!(t.nodes[par].children.contains(&child), "it {it}");
                depth += 1;
                child = par;
                cur = t.nodes[par].parent;
            }
            assert_eq!(depth, inside, "it {it}: node {i}");
            assert_eq!(t.nodes[i].is_hole, depth % 2 == 1, "it {it}: node {i}");
        }
        let ra2: Vec<TaggedRing> = a.iter().rev().cloned().collect();
        let rb2: Vec<TaggedRing> = b.iter().rev().cloned().collect();
        assert_eq!(
            run(&ra2, &rb2).to_tagged_polygons(),
            t.to_tagged_polygons(),
            "it {it}"
        );
    }
}
