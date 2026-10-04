use polyclip::*;
use proptest::prelude::*;

fn rings(range: i64, max_rings: usize, max_pts: usize) -> impl Strategy<Value = Vec<Ring>> {
    prop::collection::vec(
        prop::collection::vec((-range..=range, -range..=range), 3..max_pts)
            .prop_map(|v| v.into_iter().map(Point::from).collect::<Ring>()),
        1..max_rings,
    )
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

proptest! {
    #[test]
    fn offset_valid(a in rings(10_000, 3, 9), d in -3000i64..3000, j in join(), s in side()) {
        let r = offset(&a, d, j, ArcTol::new(20, s)).unwrap();
        prop_assert_eq!(check_canonical(&r, true), Ok(()));
        let norm = union_all(&a, FillRule::NonZero).unwrap();
        let an = area2(&norm);
        let ar = area2(&r);
        if d > 0 { prop_assert!(ar >= an); }
        if d < 0 { prop_assert!(ar <= an); }
    }

    #[test]
    fn canonical_fast_path_matches(a in rings(10_000, 3, 9), d in -3000i64..3000, j in join()) {
        // Canonical input skips normalization; rotating a ring start defeats the check
        // and must not change the result.
        let tol = ArcTol::new(20, Side::Outside);
        let norm = union_all(&a, FillRule::NonZero).unwrap();
        let fast = offset(&norm, d, j, tol).unwrap();
        let mut rot = norm.clone();
        if let Some(p) = rot.first_mut() { p.outer.rotate_left(1); }
        prop_assert_eq!(&fast, &offset(&rot, d, j, tol).unwrap());
        let tree = Boolean::new().subject(&norm, FillRule::NonZero).execute_tree().unwrap();
        prop_assert_eq!(&fast, &offset(&tree, d, j, tol).unwrap());
        if d > 0 {
            prop_assert_eq!(opening(&norm, d, tol).unwrap(), offset(&offset(&rot, -d, Join::Round, tol).unwrap(), d, Join::Round, tol).unwrap());
            prop_assert_eq!(closing(&norm, d, tol).unwrap(), offset(&offset(&rot, d, Join::Round, tol).unwrap(), -d, Join::Round, tol).unwrap());
        }
    }

    #[test]
    fn grow_then_shrink_contains(a in rings(10_000, 3, 8), d in 1i64..2000) {
        // Closing (grow then shrink with round joins) contains the original up to
        // tolerance and rounding. Rounding moves edges by < 1 unit, which near a very sharp
        // tip can retreat the tip a long way, so bound the *area* that is lost.
        let tol = ArcTol::new(5, Side::Outside);
        let norm = union_all(&a, FillRule::NonZero).unwrap();
        let g = offset(&norm, d, Join::Round, tol).unwrap();
        let back = offset(&g, -d, Join::Round, tol).unwrap();
        let lost = boolean(Op::Difference, &norm, &back, FillRule::NonZero).unwrap();
        let per: f64 = norm.iter().flat_map(|p| p.rings()).map(|r| r.edges().map(|(p, q)| {
            (((q.x - p.x) as f64).powi(2) + ((q.y - p.y) as f64).powi(2)).sqrt()
        }).sum::<f64>()).sum();
        prop_assert!((area2(&lost) as f64) / 2.0 <= 2.0 * per + 100.0, "lost {} per {}", area2(&lost), per);
    }

    #[test]
    fn path_offset_valid(paths in prop::collection::vec(prop::collection::vec((-5000i64..5000, -5000i64..5000), 1..6), 1..4),
                         d in 0i64..800, j in join(),
                         cap in prop_oneof![Just(EndCap::Round), Just(EndCap::Square), Just(EndCap::Butt), Just(EndCap::Joined)]) {
        let paths: Vec<Path> = paths.into_iter().map(|v| v.into_iter().map(Point::from).collect()).collect();
        let r = offset_paths(&paths, d, j, cap, ArcTol::new(10, Side::Nearest)).unwrap();
        prop_assert_eq!(check_canonical(&r, true), Ok(()));
        // Every path vertex is covered (for non-butt caps and d > 1).
        if d > 2 && cap != EndCap::Butt {
            for p in &paths {
                for &v in p.iter() {
                    prop_assert_ne!(locate(&r, v), Location::Outside);
                }
            }
        }
    }
}
