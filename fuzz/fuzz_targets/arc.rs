//! Arc and circle approximation (`Circle::to_ring`, `Shape::to_polygon`,
//! `Shape::to_tagged`) with every side and tolerance, degenerate arcs (collinear
//! three-point arcs, zero radius, full circles) and extreme centres / radii.
//!
//! Checks: no panic; out-of-range centres are rejected; results have in-range vertices;
//! a circle of radius >= 16 is a valid, canonically oriented ring whose approximation
//! lies on the requested side; tagged and untagged shape conversions agree and carry one
//! tag per edge; the approximated shape normalizes to a canonical polygon set.
#![no_main]

use libfuzzer_sys::fuzz_target;
use polyclip::*;
use polyclip_fuzz::{Gen, all_in_range};

fn span(pts: &[Point]) -> i64 {
    pts.iter()
        .map(|p| p.x.unsigned_abs().max(p.y.unsigned_abs()))
        .max()
        .unwrap_or(0)
        .min(i64::MAX as u64 / 4) as i64
        * 2
}

fn contour(g: &mut Gen<'_>, pts: &mut Vec<Point>) -> Contour {
    let n = g.int(0, 6);
    (0..n)
        .map(|_| {
            let e = g.point();
            pts.push(e);
            match g.int(0, 2) {
                0 => Curve::Line(e),
                1 => {
                    let mid = g.point();
                    pts.push(mid);
                    Curve::Arc { mid, end: e }
                }
                _ => {
                    let center = g.point();
                    pts.push(center);
                    Curve::CenterArc {
                        center,
                        end: e,
                        ccw: g.bool(),
                    }
                }
            }
        })
        .collect()
}

fuzz_target!(|data: &[u8]| {
    let mut g = Gen::new(data);
    if g.bool() {
        // Circle.
        let c = g.point();
        let r = match g.int(0, 3) {
            0 => g.int(-2, 3),
            1 => g.int(1, 64),
            2 => g.int(1, g.extent()),
            _ => g.int(MAX_COORD - 4, MAX_COORD + 4),
        };
        let tol = g.arc_tol_for(r);
        polyclip_fuzz::dump("(circle, tol)", &(Circle::new(c, r), tol));
        let Ok(ring) = Circle::new(c, r).to_ring(tol) else {
            return;
        };
        assert!(c.in_range() && r > 0 && tol.tolerance >= 1);
        assert!(all_in_range(ring.iter()));
        if (16..=(1 << 30)).contains(&r) {
            let poly = [Polygon::new(ring.clone(), vec![])];
            if let Err(e) = check_canonical(&poly, false) {
                panic!("circle r={r} c={c:?} {tol:?} not canonical: {e:?}");
            }
            let probe = |dx: i64| Point::new(c.x + dx, c.y);
            match tol.side {
                Side::Outside => assert_ne!(
                    locate(&ring, probe(r - 1)),
                    Location::Outside,
                    "Outside approximation misses the circle"
                ),
                Side::Inside => assert_eq!(
                    locate(&ring, probe(r + 1)),
                    Location::Outside,
                    "Inside approximation exceeds the circle"
                ),
                Side::Nearest => {}
            }
        }
    } else {
        // Shape.
        let mut pts = Vec::new();
        let outer = contour(&mut g, &mut pts);
        let nh = g.int(0, 2);
        let holes: Vec<Contour> = (0..nh).map(|_| contour(&mut g, &mut pts)).collect();
        let shape = Shape::new(outer, holes);
        let tol = g.arc_tol_for(span(&pts));
        polyclip_fuzz::dump("(shape, tol)", &(&shape, tol));
        let poly = match shape.to_polygon(tol) {
            Ok(p) => p,
            Err(_) => return,
        };
        assert!(all_in_range(pts.iter()), "out-of-range shape accepted");
        assert!(all_in_range(poly.rings().flat_map(|r| r.iter())));
        let tagged = shape
            .to_tagged(tol, &|i, j| ((i as u64) << 32) | j as u64)
            .expect("to_tagged after to_polygon succeeded");
        assert_eq!(tagged.len(), 1 + shape.holes.len());
        for (t, r) in tagged.iter().zip(poly.rings()) {
            assert_eq!(t.points, r.0, "tagged and untagged rings differ");
            if t.points.len() > 1 {
                assert_eq!(t.tags.len(), t.points.len(), "one tag per edge");
            }
            for &tag in &t.tags {
                let (ci, j) = ((tag >> 32) as usize, (tag & 0xffff_ffff) as usize);
                assert!(ci <= shape.holes.len(), "tag names a missing contour");
                let len = if ci == 0 {
                    shape.contour.len()
                } else {
                    shape.holes[ci - 1].len()
                };
                assert!(j < len, "tag names a missing element");
            }
        }
        let norm = union_all(&poly, FillRule::NonZero).expect("normalize shape");
        if let Err(e) = check_canonical(&norm, true) {
            panic!("normalized shape not canonical: {e:?}");
        }
    }
});
