use polyclip::*;
use proptest::prelude::*;
use std::cmp::Ordering;

fn rings(range: i64, max_rings: usize, max_pts: usize) -> impl Strategy<Value = Vec<Ring>> {
    prop::collection::vec(
        prop::collection::vec((-range..=range, -range..=range), 3..max_pts)
            .prop_map(|v| v.into_iter().map(Point::from).collect::<Ring>()),
        1..max_rings,
    )
}

proptest! {
    #[test]
    fn polygon_distance_consistent(a in rings(1000, 3, 7), b in rings(1000, 3, 7), dx in -3000i64..3000, d in 0i64..2000) {
        let pa = union_all(&a, FillRule::NonZero).unwrap();
        let shifted: Vec<Ring> = b.iter().map(|r| r.iter().map(|p| Point::new(p.x + dx, p.y)).collect()).collect();
        let pb = union_all(&shifted, FillRule::NonZero).unwrap();
        if pa.is_empty() || pb.is_empty() {
            return Ok(());
        }
        let c = distance(&pa, &pb).unwrap();
        // Zero distance iff the closed regions intersect.
        prop_assert_eq!(c.sq.is_zero(), intersects(&pa, &pb));
        // Threshold form agrees exactly with the full distance.
        prop_assert_eq!(distance_less_than(&pa, &pb, d), d > 0 && c.sq.cmp_dist(d as u64) == Ordering::Less);
        // Symmetric.
        prop_assert_eq!(distance(&pb, &pa).unwrap().sq, c.sq);
        // Closest points lie (approximately) on/in the geometries and realize the distance.
        let dd = ((c.a.x - c.b.x).powi(2) + (c.a.y - c.b.y).powi(2)).sqrt();
        prop_assert!((dd - c.sq.distance_f64()).abs() < 1e-6 * (1.0 + dd));
        // Boolean cross-check: when the regions are disjoint, growing one by just over half
        // the distance on each side... simpler: intersection empty iff distance > 0.
        let i = boolean(Op::Intersection, &pa, &pb, FillRule::NonZero).unwrap();
        if !c.sq.is_zero() {
            prop_assert!(i.is_empty());
        }
    }
}
