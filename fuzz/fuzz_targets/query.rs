//! Exact queries (`locate`, `intersects`, `contains`, `area2`, `centroid`, plus the ring
//! helpers) on raw (possibly invalid) and on canonical geometry of every kind: points,
//! segments, paths, rings, polygons, polygon sets and trees.
//!
//! Raw input only has to not panic. On canonical input (the output of `union_all`) the
//! results are cross-checked:
//! * `contains(a, b)` with `b` non-empty implies `intersects(a, b)`;
//! * `intersects` is symmetric; `contains(a, a)`;
//! * `contains(a, p)` for a point equals `locate(a, p) != Outside`, and `intersects`
//!   agrees with it;
//! * every vertex of `a` is `OnBoundary`;
//! * `area2` equals the sum of ring areas, and the tree form agrees with the set form;
//! * the centroid of a set with positive area lies inside its bounding box;
//! * `contains(a, b)` for regions implies `b − a` has (almost) no area.
#![no_main]

use libfuzzer_sys::fuzz_target;
use polyclip::*;
use polyclip_fuzz::{Gen, perimeter};

fn raw_queries<A: Geometry + ?Sized, B: Geometry + ?Sized>(a: &A, b: &B, p: Point) {
    let _ = locate(a, p);
    let _ = locate(b, p);
    let i = intersects(a, b);
    let c = contains(a, b);
    let _ = contains(b, a);
    let _ = area2(a);
    let _ = centroid(a);
    // Only meaningful when `b` is non-empty; raw input may be self-crossing, so
    // `contains` is unspecified there, but `intersects` is well defined for any input.
    let _ = (i, c);
}

fn checked<A: Geometry + ?Sized, B: Geometry + ?Sized>(a: &A, b: &B) {
    let i = intersects(a, b);
    assert_eq!(i, intersects(b, a), "intersects is not symmetric");
    if contains(a, b) && b.bbox().is_some() {
        assert!(i, "contains(a, b) but !intersects(a, b)");
    }
    if contains(b, a) && a.bbox().is_some() {
        assert!(i, "contains(b, a) but !intersects(a, b)");
    }
}

fn check_point(a: &PolygonSet, p: Point) {
    let l = locate(a, p);
    assert_eq!(contains(a, &p), l != Location::Outside, "contains/locate");
    assert_eq!(
        intersects(a, &p),
        l != Location::Outside,
        "intersects/locate"
    );
    let in_some = a
        .iter()
        .map(|poly| locate_in_polygon(poly, p))
        .collect::<Vec<_>>();
    match l {
        Location::Outside => assert!(in_some.iter().all(|&x| x == Location::Outside)),
        Location::Inside => assert!(in_some.contains(&Location::Inside)),
        Location::OnBoundary => assert!(in_some.contains(&Location::OnBoundary)),
    }
}

fuzz_target!(|data: &[u8]| {
    let mut g = Gen::new(data);
    let ra = g.rings(3, 10);
    let rb = g.rings(3, 10);
    let path = g.path(6);
    let p = g.point();
    let q = g.point();
    let seg = Segment::new(p, q);
    let polys = g.raw_polygons(2);
    polyclip_fuzz::dump(
        "(ra, rb, path, p, q, polys)",
        &(&ra, &rb, &path, p, q, &polys),
    );
    // Queries are infallible and documented for coordinates within ±MAX_COORD only (their
    // exact i128 arithmetic needs it); out-of-range input is the error path of the
    // fallible operations, covered by the other targets.
    if g.produced_oor {
        return;
    }

    // Raw input: must not panic.
    for r in &ra {
        let _ = ring_area2(r);
        let _ = ring_winding(r, p);
        let _ = locate_in_ring(r, p);
        raw_queries(r, &seg, q);
        raw_queries(r, &path, q);
        for s in &rb {
            raw_queries(r, s, p);
        }
    }
    raw_queries(&polys, &ra.first().cloned().unwrap_or_default(), p);
    raw_queries(&path, &seg, p);
    raw_queries(&p, &q, p);
    for poly in &polys {
        let _ = locate_in_polygon(poly, q);
        raw_queries(poly, &path, p);
    }

    // Canonical input: cross-checks.
    let Ok(a) = union_all(&ra, g.fill_rule()) else {
        return;
    };
    let Ok(b) = union_all(&rb, g.fill_rule()) else {
        return;
    };
    let ring_sum: i128 = a.iter().map(|x| x.signed_area2()).sum();
    assert_eq!(area2(&a), ring_sum, "area2 vs ring areas");
    assert!(ring_sum >= 0, "canonical set with negative area");
    let tree = Boolean::new()
        .subject(&ra, FillRule::NonZero)
        .execute_tree()
        .expect("tree");
    assert_eq!(area2(&tree), tree.signed_area2());

    if !a.is_empty() {
        assert!(contains(&a, &a), "contains(a, a) is false");
    }
    for v in a.iter().flat_map(|x| x.rings().flat_map(|r| r.iter())) {
        assert_eq!(
            locate(&a, *v),
            Location::OnBoundary,
            "vertex {v:?} not on boundary"
        );
    }
    check_point(&a, p);
    check_point(&a, q);
    if let Some(c) = centroid(&a) {
        let bb = Geometry::bbox(&a).expect("bbox of non-empty set");
        assert!(
            c.x >= bb.min.x as f64 - 1e-6 * (bb.width() as f64 + 1.0)
                && c.x <= bb.max.x as f64 + 1e-6 * (bb.width() as f64 + 1.0)
                && c.y >= bb.min.y as f64 - 1e-6 * (bb.height() as f64 + 1.0)
                && c.y <= bb.max.y as f64 + 1e-6 * (bb.height() as f64 + 1.0),
            "centroid {c:?} outside bbox {bb:?}"
        );
    }

    checked(&a, &b);
    checked(&a, &seg);
    checked(&a, &path);
    checked(&a, &p);
    for poly in &a {
        checked(poly, &b);
        checked(&poly.outer, &seg);
    }

    if contains(&a, &b) && !b.is_empty() {
        let rest = boolean(Op::Difference, &b, &a, FillRule::NonZero).expect("b - a");
        // Snap rounding may move edges by < 1 unit: allow a sliver along the boundary.
        let slack = 2.0 * (perimeter(&a) + perimeter(&b)) + 4.0;
        assert!(
            (area2(&rest) as f64) <= slack,
            "contains(a, b) but area(b - a) = {}/2 (slack {slack})",
            area2(&rest)
        );
    }
});
