use polyclip::*;
use proptest::prelude::*;

fn rings(range: i64, max_rings: usize, max_pts: usize) -> impl Strategy<Value = Vec<Ring>> {
    prop::collection::vec(
        prop::collection::vec((-range..=range, -range..=range), 3..max_pts)
            .prop_map(|v| v.into_iter().map(Point::from).collect::<Ring>()),
        1..max_rings,
    )
}

proptest! {
    #[test]
    fn fracture_exact(outer in rings(100, 3, 8), cuts in rings(100, 8, 6), small in any::<bool>()) {
        // Polygons with many holes (and holes touching things) from a difference.
        let scale = if small { 1 } else { 1000 };
        let sc = |v: &Vec<Ring>| -> Vec<Ring> { v.iter().map(|r| r.iter().map(|p| Point::new(p.x * scale, p.y * scale)).collect()).collect() };
        let polys = boolean(Op::Difference, &sc(&outer), &sc(&cuts), FillRule::EvenOdd).unwrap();
        for p in &polys {
            let f = fracture(p).unwrap();
            prop_assert_eq!(f.signed_area2(), p.signed_area2());
            let n = f.len();
            for i in 0..n {
                for j in i + 1..n {
                    let (a, b) = (f[i], f[(i + 1) % n]);
                    let (c, d) = (f[j], f[(j + 1) % n]);
                    prop_assert!(!predicates::segments_cross_properly(a, b, c, d));
                }
                prop_assert_ne!(f[i], f[(i + 1) % n]);
            }
            // Same region: the symmetric difference is empty (exact, no crossings involved).
            let x = boolean(Op::Xor, &f, p, FillRule::NonZero).unwrap();
            prop_assert!(x.is_empty(), "{:?}", x);
        }
    }
}
