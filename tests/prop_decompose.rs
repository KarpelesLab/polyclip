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
    fn trapezoids_cover_exactly(a in rings(50, 4, 8)) {
        let polys = union_all(&a, FillRule::NonZero).unwrap();
        let t = trapezoids(&polys).unwrap();
        let area: f64 = t.iter().map(|t| t.area()).sum();
        let exact = area2(&polys) as f64 / 2.0;
        prop_assert!((area - exact).abs() <= 1e-6 * (1.0 + exact.abs()), "{} vs {}", area, exact);
        for z in &t {
            prop_assert!(z.x1 > z.x0);
            let c = z.corners();
            // Positive height everywhere and the centre lies in the region (doubled coords).
            prop_assert!(c[3].y >= c[0].y && c[2].y >= c[1].y);
            let cx = (c[0].x + c[1].x + c[2].x + c[3].x) / 4.0;
            let cy = (c[0].y + c[1].y + c[2].y + c[3].y) / 4.0;
            // Sample a point slightly inside: scale by 4 to keep integer exactness approximate.
            let p = Point::new((cx * 4.0).round() as i64, (cy * 4.0).round() as i64);
            let scaled: PolygonSet = polys.iter().map(|q| Polygon {
                outer: q.outer.iter().map(|v| Point::new(4 * v.x, 4 * v.y)).collect(),
                holes: q.holes.iter().map(|h| h.iter().map(|v| Point::new(4 * v.x, 4 * v.y)).collect()).collect(),
            }).collect();
            prop_assert_ne!(locate(&scaled, p), Location::Outside);
        }
        // Disjoint interiors: overlapping x-ranges must have separated y-ranges at the
        // middle of the overlap.
        for i in 0..t.len() {
            for j in i + 1..t.len() {
                let (a, b) = (&t[i], &t[j]);
                let lo = a.x0.max(b.x0);
                let hi = a.x1.min(b.x1);
                if lo >= hi { continue; }
                let xm = (lo + hi) as f64 / 2.0;
                let y = |e: (Point, Point)| e.0.y as f64 + (xm - e.0.x as f64) * (e.1.y - e.0.y) as f64 / (e.1.x - e.0.x) as f64;
                let (ab, at, bb, bt) = (y(a.bottom), y(a.top), y(b.bottom), y(b.top));
                prop_assert!(at <= bb + 1e-9 || bt <= ab + 1e-9, "{:?} overlaps {:?}", a, b);
            }
        }
    }
}
