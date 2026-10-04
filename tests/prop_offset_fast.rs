//! Offsets through the fast paths (simple raw curves passed through without noding,
//! parallel raw curves) must be identical to the reference pipeline, for every join, cap,
//! tolerance side, delta sign, tags and holes.

use polyclip::*;
use proptest::prelude::*;

/// Runs `f` with the reference offset pipeline (and one-piece booleans), then with the
/// default one, and returns both results.
fn both<T>(f: impl Fn() -> T) -> (T, T) {
    set_offset_reference(true);
    set_always_monolithic(true);
    let a = f();
    set_offset_reference(false);
    set_always_monolithic(false);
    (a, f())
}

fn join() -> impl Strategy<Value = Join> {
    prop_oneof![
        Just(Join::Round),
        Just(Join::Bevel),
        Just(Join::Square),
        (1.0f64..4.0).prop_map(|limit| Join::Miter { limit })
    ]
}

fn side() -> impl Strategy<Value = Side> {
    prop_oneof![Just(Side::Outside), Just(Side::Inside), Just(Side::Nearest)]
}

fn cap() -> impl Strategy<Value = EndCap> {
    prop_oneof![
        Just(EndCap::Round),
        Just(EndCap::Square),
        Just(EndCap::Butt),
        Just(EndCap::Joined)
    ]
}

/// A ring: random vertices (often self-intersecting) or a convex polygon (points on a
/// circle, possibly reversed), whose raw offset curves are often simple.
fn ring(range: i64) -> impl Strategy<Value = Vec<Point>> {
    prop_oneof![
        prop::collection::vec((-range..=range, -range..=range), 3..10)
            .prop_map(|v| v.into_iter().map(Point::from).collect()),
        (3usize..40, 1i64..=range, any::<bool>(), -range..=range).prop_map(|(n, r, rev, cx)| {
            let mut v: Vec<Point> = (0..n)
                .map(|k| {
                    let a = k as f64 / n as f64 * core::f64::consts::TAU;
                    Point::new(
                        cx + (r as f64 * a.cos()).round() as i64,
                        (r as f64 * a.sin()).round() as i64,
                    )
                })
                .collect();
            if rev {
                v.reverse();
            }
            v
        }),
    ]
}

fn tagged(range: i64) -> impl Strategy<Value = Vec<TaggedRing>> {
    prop::collection::vec(
        (ring(range), any::<u64>()).prop_map(|(points, seed)| {
            let tags = (0..points.len() as u64)
                .map(|i| (seed >> (i % 60)) & 3)
                .collect();
            TaggedRing { points, tags }
        }),
        1..4,
    )
}

fn plain(rs: &[TaggedRing]) -> Vec<Ring> {
    rs.iter().map(|r| Ring(r.points.clone())).collect()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(1500))]

    #[test]
    fn offset_matches_reference(rs in tagged(5000), d in -2000i64..2000, j in join(), s in side(),
                                t in 1i64..40) {
        let tol = ArcTol::new(t, s);
        let p = plain(&rs);
        let (a, b) = both(|| offset(&p, d, j, tol));
        prop_assert_eq!(a, b);
        // A single ring (the single-curve fast path).
        let (a, b) = both(|| offset(&p[0], d, j, tol));
        prop_assert_eq!(a, b);
        let (a, b) = both(|| offset_tagged(&rs, d, j, tol, 99));
        prop_assert_eq!(a, b);
        let (a, b) = both(|| offset_tagged(&rs[0], d, j, tol, 99));
        prop_assert_eq!(a, b);
    }

    #[test]
    fn canonical_and_holes_match_reference(rs in tagged(5000), d in -2000i64..2000, j in join(),
                                           s in side()) {
        let tol = ArcTol::new(7, s);
        // Canonical input with holes (frame minus the rings), offset directly.
        let frame = Ring::from([(-6000, -6000), (6000, -6000), (6000, 6000), (-6000, 6000)]);
        let holed = boolean(Op::Difference, &frame, &plain(&rs), FillRule::NonZero).unwrap();
        let (a, b) = both(|| offset(&holed, d, j, tol));
        prop_assert_eq!(a, b);
        let tree = Boolean::new().subject(&rs, FillRule::NonZero).execute_tree().unwrap();
        let (a, b) = both(|| offset_tagged(&tree, d, j, tol, 3));
        prop_assert_eq!(a, b);
        let (a, b) = both(|| opening(&holed, d.abs(), tol));
        prop_assert_eq!(a, b);
        let (a, b) = both(|| closing(&holed, d.abs(), tol));
        prop_assert_eq!(a, b);
    }

    #[test]
    fn paths_match_reference(rs in tagged(5000), d in 0i64..1500, j in join(), c in cap(),
                             s in side()) {
        let tol = ArcTol::new(5, s);
        let paths: Vec<TaggedPath> = rs
            .iter()
            .map(|r| TaggedPath { points: r.points.clone(), tags: r.tags.clone() })
            .collect();
        let (a, b) = both(|| offset_paths_tagged(&paths, d, j, c, tol, 8));
        prop_assert_eq!(a, b);
        let (a, b) = both(|| offset_paths_tagged(&paths[0], d, j, c, tol, 8));
        prop_assert_eq!(a, b);
        // A two-point track (the most common stroke).
        let track = Path(paths[0].points[..2].to_vec());
        let (a, b) = both(|| offset_paths(&track, d, j, c, tol));
        prop_assert_eq!(a, b);
    }

    #[test]
    fn small_coordinates_match_reference(rs in tagged(12), d in -8i64..8, j in join(),
                                         s in side(), c in cap()) {
        // Dense integer grids: rounding, collinear vertices and touching curves.
        let tol = ArcTol::new(1, s);
        let p = plain(&rs);
        let (a, b) = both(|| offset(&p[0], d, j, tol));
        prop_assert_eq!(a, b);
        let (a, b) = both(|| offset_tagged(&rs, d, j, tol, 1));
        prop_assert_eq!(a, b);
        let path = Path(p[0].0.clone());
        let (a, b) = both(|| offset_paths(&path, d.abs(), j, c, tol));
        prop_assert_eq!(a, b);
    }
}

/// Large inputs, so that raw curves are computed in parallel chunks with `rayon`.
#[test]
fn large_inputs_match_reference() {
    let tol = ArcTol::new(1000, Side::Outside);
    let circle = |cx: i64, r: f64, n: usize| -> Ring {
        (0..n)
            .map(|k| {
                let a = k as f64 / n as f64 * core::f64::consts::TAU;
                Point::new(
                    cx + (r * a.cos()).round() as i64,
                    (r * a.sin()).round() as i64,
                )
            })
            .collect()
    };
    let big = circle(0, 50_000_000.0, 10_000);
    for d in [100_000, -100_000, 2] {
        for j in [Join::Round, Join::Miter { limit: 2.0 }, Join::Square] {
            let (a, b) = both(|| offset(&big, d, j, tol).unwrap());
            assert_eq!(a, b);
        }
    }
    let many: Vec<Ring> = (0..3000)
        .map(|k| circle(k * 3_000_000, 1_000_000.0 + (k % 7) as f64 * 1e5, 24))
        .collect();
    let set = union_all(&many, FillRule::NonZero).unwrap();
    for d in [400_000, -300_000] {
        let (a, b) = both(|| offset(&set, d, Join::Round, tol).unwrap());
        assert_eq!(a, b);
    }
}
