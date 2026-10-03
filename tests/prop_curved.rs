//! Property tests for arc-preserving booleans on curved shapes (`curved_boolean`).
//!
//! Random operands made of circles, rounded rectangles, three-point arc segments,
//! triangles and annuli (often snapped to a coarse grid so that tangencies, coincident arcs
//! and shared vertices are frequent). For every case:
//!
//! * the output re-approximated finely (tolerance 1) and at the operation tolerance
//!   (`Side::Nearest`, the certified approximation) is a valid polygon set;
//! * its area matches the polygon boolean of fine approximations within the documented
//!   band bound (perimeter x 2 (t + 2), plus approximation slack);
//! * every output arc lies on an input circle (same centre, end points within rounding of
//!   one of the input radii for that centre);
//! * the result agrees with the exact result at sample points farther than `t + 6` from
//!   every input boundary, and the side guarantee holds at points farther than 6: with
//!   `Side::Inside` no point outside the exact result is inside the output (e.g. a zone
//!   minus obstacles never intrudes into a true obstacle), with `Side::Outside` no point
//!   of the exact result is missing.

use polyclip::*;
use proptest::prelude::*;
use std::f64::consts::PI;

fn p(x: i64, y: i64) -> Point {
    Point::new(x, y)
}

#[derive(Clone, Debug)]
enum Gen {
    Circle(i64, i64, i64),
    Rounded(i64, i64, i64, i64, i64),
    Dee(i64, i64, i64, i64, i64),
    Tri(i64, i64, i64, i64, i64, i64),
    Annulus(i64, i64, i64, i64),
}

fn circle_contour(c: Point, r: i64, ccw: bool) -> Contour {
    vec![Curve::CenterArc {
        center: c,
        end: p(c.x + r, c.y),
        ccw,
    }]
}

fn build(g: &Gen) -> Shape {
    match *g {
        Gen::Circle(x, y, r) => Shape::new(circle_contour(p(x, y), r, true), vec![]),
        Gen::Annulus(x, y, r, k) => {
            let inner = (r * k / 10).max(1);
            Shape::new(
                circle_contour(p(x, y), r, true),
                vec![circle_contour(p(x + (r - inner) / 3, y), inner, false)],
            )
        }
        Gen::Rounded(x0, y0, w, h, rr) => {
            let (x1, y1) = (x0 + w, y0 + h);
            let r = rr.min(w / 2).min(h / 2).max(1);
            Shape::new(
                vec![
                    Curve::Line(p(x1 - r, y0)),
                    Curve::CenterArc {
                        center: p(x1 - r, y0 + r),
                        end: p(x1, y0 + r),
                        ccw: true,
                    },
                    Curve::Line(p(x1, y1 - r)),
                    Curve::CenterArc {
                        center: p(x1 - r, y1 - r),
                        end: p(x1 - r, y1),
                        ccw: true,
                    },
                    Curve::Line(p(x0 + r, y1)),
                    Curve::CenterArc {
                        center: p(x0 + r, y1 - r),
                        end: p(x0, y1 - r),
                        ccw: true,
                    },
                    Curve::Line(p(x0, y0 + r)),
                    Curve::CenterArc {
                        center: p(x0 + r, y0 + r),
                        end: p(x0 + r, y0),
                        ccw: true,
                    },
                ],
                vec![],
            )
        }
        // A circular segment: chord from (x, y) to (x + w, y) and an arc bulging by `b`
        // (three-point form, usually with a non-integral centre).
        Gen::Dee(x, y, w, b, flip) => {
            let b = if flip % 2 == 0 { b } else { -b };
            Shape::new(
                vec![
                    Curve::Line(p(x + w, y)),
                    Curve::Arc {
                        mid: p(x + w / 2, y + b),
                        end: p(x, y),
                    },
                ],
                vec![],
            )
        }
        Gen::Tri(a, b, c, d, e, f) => Shape::new(
            vec![
                Curve::Line(p(c, d)),
                Curve::Line(p(e, f)),
                Curve::Line(p(a, b)),
            ],
            vec![],
        ),
    }
}

fn gen_shape() -> impl Strategy<Value = Gen> {
    // Coarse grid (q = 100 000) makes tangencies and coincidences likely; fine grid is
    // generic.
    (prop_oneof![Just(1i64), Just(1000), Just(100_000)], 0..5u8).prop_flat_map(|(q, kind)| {
        let c = -10..=10i64;
        let s = 1..=8i64;
        match kind {
            0 => (c.clone(), c.clone(), s.clone())
                .prop_map(move |(x, y, r)| {
                    let q2 = q.max(1000);
                    Gen::Circle(
                        x * 100_000 / q * q,
                        y * 100_000 / q * q,
                        r * 60_000 / q2 * q2,
                    )
                })
                .boxed(),
            1 => (c.clone(), c.clone(), s.clone(), s.clone(), 1..=6i64)
                .prop_map(move |(x, y, w, h, r)| {
                    Gen::Rounded(
                        x * 100_000 / q * q,
                        y * 100_000 / q * q,
                        w * 100_000,
                        h * 100_000,
                        r * 50_000,
                    )
                })
                .boxed(),
            2 => (c.clone(), c.clone(), s.clone(), 1..=6i64, 0..2i64)
                .prop_map(move |(x, y, w, b, f)| {
                    Gen::Dee(
                        x * 100_000 + q % 7,
                        y * 100_000,
                        w * 100_000,
                        b * 40_000 + 1,
                        f,
                    )
                })
                .boxed(),
            3 => (
                c.clone(),
                c.clone(),
                c.clone(),
                c.clone(),
                c.clone(),
                c.clone(),
            )
                .prop_map(move |(a, b, cc, d, e, f)| {
                    let k = 100_000 / q * q;
                    Gen::Tri(a * k + 3, b * k, cc * k, d * k + 7, e * k, f * k)
                })
                .prop_filter("non-degenerate", |g| match *g {
                    Gen::Tri(a, b, c, d, e, f) => (c - a) * (f - b) - (d - b) * (e - a) != 0,
                    _ => true,
                })
                .boxed(),
            _ => (c.clone(), c, s, 1..=8i64)
                .prop_map(move |(x, y, r, k)| {
                    Gen::Annulus(x * 100_000 / q * q, y * 100_000, r * 60_000, k)
                })
                .boxed(),
        }
    })
}

fn op() -> impl Strategy<Value = Op> {
    prop_oneof![
        Just(Op::Union),
        Just(Op::Intersection),
        Just(Op::Difference),
        Just(Op::Xor)
    ]
}

fn side() -> impl Strategy<Value = Side> {
    prop_oneof![Just(Side::Inside), Just(Side::Outside), Just(Side::Nearest)]
}

fn approx_all(s: &[Shape], tol: ArcTol) -> Vec<Polygon> {
    s.iter().map(|x| x.to_polygon(tol).unwrap()).collect()
}

fn area(ps: &[Polygon]) -> f64 {
    ps.iter().map(|p| p.signed_area2() as f64 / 2.0).sum()
}

fn seg_dist(q: (f64, f64), a: Point, b: Point) -> f64 {
    let (ax, ay, bx, by) = (a.x as f64, a.y as f64, b.x as f64, b.y as f64);
    let (dx, dy) = (bx - ax, by - ay);
    let l2 = dx * dx + dy * dy;
    let t = if l2 > 0.0 {
        (((q.0 - ax) * dx + (q.1 - ay) * dy) / l2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    ((q.0 - ax - t * dx).powi(2) + (q.1 - ay - t * dy).powi(2)).sqrt()
}

fn boundary_dist(poly: &Polygon, q: Point) -> f64 {
    let qf = (q.x as f64, q.y as f64);
    let mut d = f64::MAX;
    for r in poly.rings() {
        for (a, b) in r.edges() {
            d = d.min(seg_dist(qf, a, b));
        }
    }
    d
}

fn perimeter(ps: &[Polygon]) -> f64 {
    ps.iter()
        .flat_map(|p| p.rings())
        .flat_map(|r| r.edges())
        .map(|(a, b)| seg_dist((a.x as f64, a.y as f64), b, b))
        .sum()
}

/// Input circles: centre -> radii (three-point arcs with their rounded circumcentre).
fn source_circles(shapes: &[Shape]) -> Vec<(Point, f64)> {
    let mut out = Vec::new();
    for s in shapes {
        for c in std::iter::once(&s.contour).chain(s.holes.iter()) {
            let mut cur = c.last().map(|e| e.end()).unwrap_or_default();
            for e in c {
                match *e {
                    Curve::CenterArc { center, .. } => {
                        let r = (((cur.x - center.x) as f64).powi(2)
                            + ((cur.y - center.y) as f64).powi(2))
                        .sqrt();
                        out.push((center, r));
                    }
                    Curve::Arc { mid, end } => {
                        let (ax, ay) = (cur.x as f64, cur.y as f64);
                        let (bx, by) = (mid.x as f64 - ax, mid.y as f64 - ay);
                        let (cx, cy) = (end.x as f64 - ax, end.y as f64 - ay);
                        let d = 2.0 * (bx * cy - by * cx);
                        if d != 0.0 {
                            let b2 = bx * bx + by * by;
                            let c2 = cx * cx + cy * cy;
                            let ux = (cy * b2 - by * c2) / d;
                            let uy = (bx * c2 - cx * b2) / d;
                            let center = p((ax + ux).round() as i64, (ay + uy).round() as i64);
                            out.push((center, (ux * ux + uy * uy).sqrt()));
                        }
                    }
                    Curve::Line(_) => {}
                }
                cur = e.end();
            }
        }
    }
    out
}

fn check_arcs(out: &[Shape], src: &[(Point, f64)]) -> std::result::Result<usize, TestCaseError> {
    let mut n = 0;
    for s in out {
        for c in std::iter::once(&s.contour).chain(s.holes.iter()) {
            let mut cur = c.last().map(|e| e.end()).unwrap_or_default();
            for e in c {
                if let Curve::CenterArc { center, end, .. } = *e {
                    n += 1;
                    let d = |q: Point| {
                        (((q.x - center.x) as f64).powi(2) + ((q.y - center.y) as f64).powi(2))
                            .sqrt()
                    };
                    let ok = src.iter().any(|&(c, r)| {
                        c == center && (d(cur) - r).abs() <= 2.0 && (d(end) - r).abs() <= 2.0
                    });
                    prop_assert!(
                        ok,
                        "arc {cur:?} -> {end:?} around {center:?} not on a source"
                    );
                }
                cur = e.end();
            }
        }
    }
    Ok(n)
}

fn inside_op(op: Op, a: bool, b: bool) -> bool {
    match op {
        Op::Union => a || b,
        Op::Intersection => a && b,
        Op::Difference => a && !b,
        Op::Xor => a != b,
    }
}

fn check_case(
    sa: &[Shape],
    sb: &[Shape],
    op: Op,
    side: Side,
    t: i64,
    samples: &[(i64, i64)],
) -> std::result::Result<(), TestCaseError> {
    let tol = ArcTol::new(t, side);
    let out = curved_boolean(op, sa, sb, FillRule::NonZero, tol).unwrap();

    // Validity, finely and at the operation tolerance.
    let fine = approx_all(&out, ArcTol::new(1, Side::Nearest));
    prop_assert!(
        validate_set(&fine).is_ok(),
        "fine: {:?}\n{out:?}",
        validate_set(&fine)
    );
    // (One-sided coarse approximations may legitimately self-intersect at sharp corners,
    // e.g. where two curves meet tangentially, so only the certified one is checked.)
    let coarse = approx_all(&out, ArcTol::new(t, Side::Nearest));
    prop_assert!(validate_set(&coarse).is_ok(), "{:?}", validate_set(&coarse));

    // Area against the polygon boolean of fine approximations.
    let fa = approx_all(sa, ArcTol::new(1, Side::Nearest));
    let fb = approx_all(sb, ArcTol::new(1, Side::Nearest));
    let reference = boolean(op, &fa, &fb, FillRule::NonZero).unwrap();
    let per = perimeter(&fa) + perimeter(&fb);
    let bound = per * 2.0 * (t as f64 + 4.0) + 100.0;
    let (ac, ar) = (area(&fine), area(&reference));
    prop_assert!(
        (ac - ar).abs() <= bound,
        "area {ac} vs {ar} (bound {bound})"
    );

    // Arcs on source circles.
    let src = source_circles(&sa.iter().chain(sb.iter()).cloned().collect::<Vec<_>>());
    check_arcs(&out, &src)?;

    // Exact membership at sample points far from every input boundary.
    let fine_in = |polys: &[Polygon], q: Point| -> (bool, f64) {
        let mut inside = false;
        let mut d = f64::MAX;
        for pg in polys {
            d = d.min(boundary_dist(pg, q));
            if locate_in_polygon(pg, q) == Location::Inside {
                inside = true;
            }
        }
        (inside, d)
    };
    let mut pts: Vec<Point> = samples.iter().map(|&(x, y)| p(x, y)).collect();
    // Also points just inside and outside every circle of the clip operand.
    for &(c, r) in &source_circles(&sa.iter().chain(sb.iter()).cloned().collect::<Vec<_>>()) {
        for k in 0..8 {
            let ang = k as f64 * PI / 4.0 + 0.3;
            for dr in [-8.0, 8.0, -(t as f64) - 8.0] {
                let rr = r + dr;
                if rr > 0.0 {
                    pts.push(p(
                        c.x + (rr * ang.cos()).round() as i64,
                        c.y + (rr * ang.sin()).round() as i64,
                    ));
                }
            }
        }
    }
    for q in pts {
        let (ina, da) = fine_in(&fa, q);
        let (inb, db) = fine_in(&fb, q);
        let d = da.min(db);
        if d <= 6.0 {
            continue;
        }
        let truth = inside_op(op, ina, inb);
        let got = fine
            .iter()
            .any(|pg| locate_in_polygon(pg, q) == Location::Inside);
        if d > t as f64 + 6.0 {
            prop_assert_eq!(got, truth, "far point {:?} (d = {})", q, d);
        }
        match side {
            Side::Inside => prop_assert!(
                !got || truth,
                "{:?} outside the result but in the output (d = {})",
                q,
                d
            ),
            Side::Outside => prop_assert!(
                got || !truth,
                "{:?} in the result but not in the output (d = {})",
                q,
                d
            ),
            Side::Nearest => {}
        }
    }
    Ok(())
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(
        std::env::var("PROPTEST_CASES").ok().and_then(|s| s.parse().ok()).unwrap_or(200)
    ))]

    #[test]
    fn curved_boolean_properties(
        a in prop::collection::vec(gen_shape(), 1..4),
        b in prop::collection::vec(gen_shape(), 0..4),
        op in op(),
        side in side(),
        t in prop_oneof![Just(10i64), Just(200), Just(2_000)],
        samples in prop::collection::vec((-1_200_000i64..1_200_000, -1_200_000i64..1_200_000), 60),
    ) {
        let sa: Vec<Shape> = a.iter().map(build).collect();
        let sb: Vec<Shape> = b.iter().map(build).collect();
        check_case(&sa, &sb, op, side, t, &samples)?;
    }

    /// Nearly tangent, nearly coincident and nearly concentric circles and arcs, offset by
    /// a few multiples of the tolerance (the hard cases for reconstruction).
    #[test]
    fn near_degenerate(
        base in gen_shape(),
        rel in prop::collection::vec((0..4u8, -40i64..40, 0..8i64, -3i64..3), 1..4),
        op in op(),
        side in side(),
        t in prop_oneof![Just(10i64), Just(200)],
        samples in prop::collection::vec((-1_200_000i64..1_200_000, -1_200_000i64..1_200_000), 40),
    ) {
        let sa = vec![build(&base)];
        // Anchor circle of the base shape (its first arc, or a default).
        let src = source_circles(&sa);
        let (c, r) = src.first().copied().unwrap_or((p(0, 0), 300_000.0));
        let r = r.round() as i64;
        let mut sb = Vec::new();
        for &(kind, d, k, s) in &rel {
            let d = d * t / 8 + s;
            let r2 = (r * (k + 2) / 8).max(10);
            let g = match kind {
                // Externally tangent (+- d).
                0 => Gen::Circle(c.x + r + r2 + d, c.y + s, r2),
                // Internally tangent.
                1 => Gen::Circle(c.x + r - r2 + d, c.y, r2),
                // Concentric, radius +- d.
                2 => Gen::Circle(c.x + s, c.y, (r + d).max(5)),
                // A chord line grazing the circle (triangle edge at y = c.y + r + d).
                _ => Gen::Tri(c.x - 3 * r, c.y + r + d, c.x + 3 * r, c.y + r + d + s, c.x, c.y + 2 * r + 50),
            };
            sb.push(build(&g));
        }
        check_case(&sa, &sb, op, side, t, &samples)?;
    }

    /// The zone-fill case: zone minus obstacles with `Side::Inside` never intrudes into an
    /// obstacle, and obstacle circles come out as full circles when isolated.
    #[test]
    fn zone_minus_obstacles_side(
        obstacles in prop::collection::vec((-20..20i64, -20..20i64, 1..6i64, 0..3i64), 1..12),
        t in prop_oneof![Just(50i64), Just(1_000)],
    ) {
        let zone = build(&Gen::Rounded(-1_000_000, -1_000_000, 2_000_000, 2_000_000, 300_000));
        let obs: Vec<Shape> = obstacles
            .iter()
            .map(|&(x, y, r, j)| build(&Gen::Circle(x * 50_000 + j * 7, y * 50_000, r * 40_000 + j)))
            .collect();
        let out = curved_boolean(
            Op::Difference,
            &[zone],
            &obs,
            FillRule::NonZero,
            ArcTol::new(t, Side::Inside),
        )
        .unwrap();
        let fine = approx_all(&out, ArcTol::new(1, Side::Outside));
        prop_assert!(validate_set(&fine).is_ok());
        for &(x, y, r, j) in &obstacles {
            let (cx, cy, rr) = ((x * 50_000 + j * 7) as f64, (y * 50_000) as f64, (r * 40_000 + j) as f64);
            for k in 0..64 {
                let ang = k as f64 * PI / 32.0;
                for depth in [5.0, 50.0, rr / 2.0, rr] {
                    let q = p(
                        (cx + (rr - depth) * ang.cos()).round() as i64,
                        (cy + (rr - depth) * ang.sin()).round() as i64,
                    );
                    let hit = fine.iter().any(|pg| locate_in_polygon(pg, q) == Location::Inside);
                    prop_assert!(!hit, "{q:?} inside obstacle ({cx}, {cy}, {rr}) is filled");
                }
            }
        }
        let src = source_circles(&obs);
        check_arcs(&out, &src.iter().copied().chain([
            (p(-700_000, -700_000), 300_000.0), (p(700_000, -700_000), 300_000.0),
            (p(700_000, 700_000), 300_000.0), (p(-700_000, 700_000), 300_000.0),
        ]).collect::<Vec<_>>())?;
    }
}
