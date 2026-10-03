//! Exact distances: `distance`, `distance_sq` and `distance_less_than` on every geometry
//! kind. Checks: symmetry, zero distance exactly when the geometries intersect, the
//! threshold form agrees with the full distance, and the closest points realize it.
#![no_main]

use libfuzzer_sys::fuzz_target;
use polyclip::*;
use polyclip_fuzz::{Gen, all_in_range, poly_points};
use std::cmp::Ordering;

fn check<A: Geometry + ?Sized, B: Geometry + ?Sized>(a: &A, b: &B, d: i64) {
    let Some(c) = distance(a, b) else {
        assert!(!distance_less_than(a, b, d));
        return;
    };
    assert_eq!(
        distance(b, a).map(|c| c.sq),
        Some(c.sq),
        "distance not symmetric"
    );
    assert_eq!(
        c.sq.is_zero(),
        intersects(a, b),
        "zero distance iff intersecting"
    );
    let lt = d > 0 && c.sq.cmp_dist(d as u64) == Ordering::Less;
    assert_eq!(
        distance_less_than(a, b, d),
        lt,
        "threshold disagrees with distance"
    );
    let dd = ((c.a.x - c.b.x).powi(2) + (c.a.y - c.b.y).powi(2)).sqrt();
    let df = c.sq.distance_f64();
    assert!(
        (dd - df).abs() <= 1e-6 * (1.0 + df) + 1e-3,
        "closest points {dd} vs {df}"
    );
}

fuzz_target!(|data: &[u8]| {
    let mut g = Gen::new(data);
    let ra = g.raw_polygons(2);
    let rb = g.raw_polygons(2);
    let path = g.path(6);
    let p = g.point();
    let d = g.int(0, 1 << 20);
    polyclip_fuzz::dump("(ra, rb, path, p, d)", &(&ra, &rb, &path, p, d));
    if !all_in_range(
        poly_points(&ra)
            .chain(poly_points(&rb))
            .chain(path.iter())
            .chain([&p]),
    ) {
        // Out-of-range input: still no panic.
        let _ = distance_less_than(&ra, &rb, d);
        return;
    }
    let (Ok(a), Ok(b)) = (
        union_all(&ra, FillRule::NonZero),
        union_all(&rb, FillRule::NonZero),
    ) else {
        return;
    };
    check(&a, &b, d);
    check(&a, &path, d);
    check(&path, &p, d);
    check(&a, &p, d);
    let seg = Segment::new(p, path.first().copied().unwrap_or(p));
    check(&seg, &b, d);
});
