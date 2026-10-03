//! Fracture, triangulation (plain and Delaunay), simplification, trapezoidal
//! decomposition, convex hull and Minkowski sums.
//!
//! Raw (possibly invalid) input must not panic. On canonical input (the output of
//! `union_all`) the results are checked exactly:
//! * fracture: one outline per polygon with the same signed area, no proper crossings;
//! * triangulation: counter-clockwise non-degenerate triangles whose areas sum to the
//!   region's area;
//! * simplification: the output is valid, keeps the polygon count, and never adds vertices;
//! * trapezoids: positive widths and total area equal to the region's (within float error);
//! * convex hull: contains every vertex;
//! * Minkowski sum with a small square: contains the polygon.
#![no_main]

use libfuzzer_sys::fuzz_target;
use polyclip::*;
use polyclip_fuzz::{Gen, all_in_range, poly_points};

fuzz_target!(|data: &[u8]| {
    let mut g = Gen::new(data);
    let raw = g.raw_polygons(3);
    let tol = g.int(0, 2000);
    polyclip_fuzz::dump("(raw, tol)", &(&raw, tol));
    // Raw input: no panics.
    for p in &raw {
        let _ = fracture(p);
        let _ = triangulate(p);
        let _ = triangulate_delaunay(p);
        let _ = simplify_polygon(p, tol);
    }
    let _ = simplify_polygons(&raw, tol);
    let _ = trapezoids(&raw);
    let _ = convex_hull(poly_points(&raw).copied());
    if !all_in_range(poly_points(&raw)) {
        return;
    }
    let Ok(set) = union_all(&raw, FillRule::EvenOdd) else {
        return;
    };
    let total = area2(&set);
    for p in &set {
        // Fracture.
        let f = fracture(p).expect("fracture of canonical polygon");
        assert_eq!(
            f.signed_area2(),
            p.signed_area2(),
            "fracture changed the area"
        );
        let n = f.len();
        for i in 0..n {
            for j in i + 1..n {
                let (a, b) = (f[i], f[(i + 1) % n]);
                let (c, d) = (f[j], f[(j + 1) % n]);
                assert!(
                    !predicates::segments_cross_properly(a, b, c, d),
                    "fracture crosses itself"
                );
            }
        }
        // Triangulation.
        for t in [
            triangulate(p).expect("triangulate"),
            triangulate_delaunay(p).expect("delaunay"),
        ] {
            let mut a2: i128 = 0;
            for tri in &t.triangles {
                let [a, b, c] = tri.map(|i| t.vertices[i as usize]);
                let o = predicates::orient(a, b, c);
                assert!(o > 0, "degenerate or clockwise triangle");
                a2 += o;
            }
            assert_eq!(a2, p.signed_area2(), "triangulation area");
        }
        // Minkowski sum with a small square contains the polygon.
        let sq = Polygon::from(Ring::from([(-2, -2), (2, -2), (2, 2), (-2, 2)]));
        if let Ok(m) = minkowski_sum(p, &sq) {
            assert!(
                contains(&m, p),
                "minkowski sum does not contain the polygon"
            );
        }
        // Convex hull.
        let h = convex_hull(p.outer.iter().copied()).expect("hull of in-range points");
        if h.len() >= 3 {
            for &v in p.outer.iter() {
                assert_ne!(locate(&h, v), Location::Outside, "hull misses a vertex");
            }
        }
    }
    // Simplification.
    let s = simplify_polygons(&set, tol);
    assert_eq!(validate_set(&s), Ok(()), "simplified set invalid");
    assert_eq!(s.len(), set.len());
    assert!(
        s.iter().map(|p| p.vertex_count()).sum::<usize>()
            <= set.iter().map(|p| p.vertex_count()).sum::<usize>()
    );
    // Trapezoids.
    let t = trapezoids(&set).expect("trapezoids");
    let a: f64 = t.iter().map(|z| z.area()).sum();
    let exact = total as f64 / 2.0;
    assert!(
        (a - exact).abs() <= 1e-6 * (1.0 + exact.abs()),
        "trapezoid area {a} vs {exact}"
    );
    assert!(t.iter().all(|z| z.x1 > z.x0));
});
