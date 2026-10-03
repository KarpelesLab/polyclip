//! End-to-end check of a PCB-like zone fill: board outline with rounded corners, pads,
//! vias and 45/90 degree tracks, clearances, the copper pour, minimum-width opening,
//! Gerber fracture, triangulation and DRC distance checks. Units are nanometers.

use polyclip::*;

const MM: i64 = 1_000_000;

fn lcg(s: &mut u64) -> u64 {
    *s = s
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    *s >> 33
}

/// A rounded rectangle pad centred at `c`.
fn rounded_rect(c: Point, w: i64, h: i64, r: i64) -> Shape {
    let (x0, y0, x1, y1) = (c.x - w / 2, c.y - h / 2, c.x + w / 2, c.y + h / 2);
    let p = Point::new;
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

struct Board {
    outline: Shape,
    pads: Vec<Shape>,
    vias: Vec<Circle>,
    tracks: Vec<Path>,
}

fn board(seed: u64) -> Board {
    let mut s = seed;
    let outline = rounded_rect(Point::new(25 * MM, 20 * MM), 50 * MM, 40 * MM, 3 * MM);
    let mut pads = Vec::new();
    // Two rows of fine-pitch pads (an IC footprint) and scattered passives.
    for i in 0..20 {
        pads.push(rounded_rect(
            Point::new(10 * MM + i * 650_000, 18 * MM),
            400_000,
            1_500_000,
            100_000,
        ));
        pads.push(rounded_rect(
            Point::new(10 * MM + i * 650_000, 24 * MM),
            400_000,
            1_500_000,
            100_000,
        ));
    }
    for _ in 0..40 {
        let c = Point::new(
            3 * MM + (lcg(&mut s) % 44_000_000) as i64,
            3 * MM + (lcg(&mut s) % 34_000_000) as i64,
        );
        pads.push(rounded_rect(c, 1_000_000, 1_200_000, 250_000));
    }
    let vias: Vec<Circle> = (0..60)
        .map(|_| {
            let c = Point::new(
                3 * MM + (lcg(&mut s) % 44_000_000) as i64,
                3 * MM + (lcg(&mut s) % 34_000_000) as i64,
            );
            Circle::new(c, 300_000)
        })
        .collect();
    // Tracks with 45 and 90 degree segments.
    let tracks: Vec<Path> = (0..30)
        .map(|_| {
            let mut p = Point::new(
                3 * MM + (lcg(&mut s) % 40_000_000) as i64,
                3 * MM + (lcg(&mut s) % 30_000_000) as i64,
            );
            let mut v = vec![p];
            for _ in 0..4 {
                let len = (lcg(&mut s) % 5_000_000) as i64 + MM;
                let dir = lcg(&mut s) % 8;
                let (dx, dy) = [
                    (1, 0),
                    (1, 1),
                    (0, 1),
                    (-1, 1),
                    (-1, 0),
                    (-1, -1),
                    (0, -1),
                    (1, -1),
                ][dir as usize];
                p = Point::new(
                    (p.x + dx * len).clamp(2 * MM, 48 * MM),
                    (p.y + dy * len).clamp(2 * MM, 38 * MM),
                );
                v.push(p);
            }
            Path(v)
        })
        .collect();
    Board {
        outline,
        pads,
        vias,
        tracks,
    }
}

#[test]
fn zone_fill_pipeline() {
    let b = board(7);
    let clearance = 200_000;
    let track_half_width = 125_000;
    let out_tol = ArcTol::new(1_000, Side::Outside);
    let in_tol = ArcTol::new(1_000, Side::Inside);

    // Obstacles: everything inflated by the clearance, approximated outward.
    let mut obstacles: Vec<Polygon> = Vec::new();
    for p in &b.pads {
        obstacles.extend(offset_shape(p, clearance, Join::Round, out_tol).unwrap());
    }
    for v in &b.vias {
        obstacles.extend(
            offset(
                &v.to_ring(out_tol).unwrap(),
                clearance,
                Join::Round,
                out_tol,
            )
            .unwrap(),
        );
    }
    obstacles.extend(
        offset_paths(
            &b.tracks,
            track_half_width + clearance,
            Join::Round,
            EndCap::Round,
            out_tol,
        )
        .unwrap(),
    );
    for p in &obstacles {
        assert_eq!(validate(p), Ok(()));
    }

    // Board outline approximated inward, shrunk by the edge clearance.
    let outline = b.outline.to_polygon(in_tol).unwrap();
    let area_outline = area2(&outline) as f64 / 2.0;
    let true_outline =
        50.0 * 40.0 * (MM * MM) as f64 - (4.0 - std::f64::consts::PI) * (3 * MM * 3 * MM) as f64;
    assert!(area_outline <= true_outline && area_outline > true_outline * 0.9999);
    let keepin = offset(&outline, -clearance, Join::Round, in_tol).unwrap();

    // The pour.
    let fill = Boolean::new()
        .subject(&keepin, FillRule::NonZero)
        .clip(&obstacles, FillRule::NonZero)
        .op(Op::Difference)
        .execute_tree()
        .unwrap();
    let fill_set = fill.to_polygon_set();
    assert_eq!(check_canonical(&fill_set, true), Ok(()));
    assert!(!fill_set.is_empty());

    // Minimum width 0.3 mm.
    let fill2 = opening(&fill_set, 150_000, in_tol).unwrap();
    assert_eq!(check_canonical(&fill2, true), Ok(()));
    assert!(area2(&fill2) <= area2(&fill_set));
    // Remove islands smaller than 1 mm^2.
    let fill3: PolygonSet = fill2
        .into_iter()
        .filter(|p| p.signed_area2() > 2 * MM as i128 * MM as i128)
        .collect();

    // DRC: the pour keeps its clearance to every obstacle source (up to rounding).
    for v in &b.vias {
        let via = v.to_ring(in_tol).unwrap();
        assert!(
            !distance_less_than(&fill3, &via, clearance - 2),
            "via clearance violated"
        );
    }
    for p in &b.pads {
        let pad = p.to_polygon(in_tol).unwrap();
        assert!(
            !distance_less_than(&fill3, &pad, clearance - 2),
            "pad clearance violated"
        );
    }
    for t in &b.tracks {
        assert!(
            !distance_less_than(&fill3, t, clearance + track_half_width - 2),
            "track clearance violated"
        );
    }
    assert!(contains(&outline, &fill3));

    // Gerber regions: fractured outlines enclose exactly the pour.
    let regions = fracture_set(&fill3).unwrap();
    let back = union_all(&regions, FillRule::NonZero).unwrap();
    assert!(
        boolean(Op::Xor, &back, &fill3, FillRule::NonZero)
            .unwrap()
            .is_empty()
    );
    let a_regions: i128 = regions.iter().map(|r| r.signed_area2()).sum();
    assert_eq!(a_regions, area2(&fill3));

    // 3D / rendering: triangulation covers the pour exactly.
    let tri = triangulate_set(&fill3).unwrap();
    let mut a_tri: i128 = 0;
    for t in &tri.triangles {
        let [a, b2, c] = t.map(|i| tri.vertices[i as usize]);
        let o = predicates::orient(a, b2, c);
        assert!(o > 0);
        a_tri += o;
    }
    assert_eq!(a_tri, area2(&fill3));

    // Arc provenance survives the offset: tag pad arcs and check the obstacle keeps them.
    let pad = &b.pads[0];
    let tagged = offset_shape_tagged(
        pad,
        clearance,
        Join::Round,
        out_tol,
        &|_, j| j as u64 + 1,
        0,
    )
    .unwrap();
    let tags: std::collections::BTreeSet<u64> = tagged.nodes[0].tags.iter().copied().collect();
    for j in 1..=8u64 {
        assert!(tags.contains(&j), "element {j} lost its tag: {tags:?}");
    }

    // Simplification keeps the pour valid.
    let simple = simplify_polygons(&fill3, 5_000);
    assert_eq!(validate_set(&simple), Ok(()));
    assert!(
        simple.iter().map(|p| p.vertex_count()).sum::<usize>()
            < fill3.iter().map(|p| p.vertex_count()).sum::<usize>()
    );
}

#[test]
fn deterministic() {
    let b = board(11);
    let tol = ArcTol::new(1_000, Side::Outside);
    let run = |reverse: bool| {
        let mut obstacles: Vec<Polygon> = Vec::new();
        for p in &b.pads {
            obstacles.extend(offset_shape(p, 200_000, Join::Round, tol).unwrap());
        }
        obstacles
            .extend(offset_paths(&b.tracks, 300_000, Join::Round, EndCap::Round, tol).unwrap());
        if reverse {
            obstacles.reverse();
            for p in obstacles.iter_mut() {
                // Also rotate every ring's start vertex.
                p.outer.rotate_left(1);
            }
        }
        let outline = b.outline.to_polygon(tol).unwrap();
        boolean(Op::Difference, &outline, &obstacles, FillRule::NonZero).unwrap()
    };
    // Input order and ring start vertices must not change the result.
    assert_eq!(run(false), run(true));
}

#[test]
fn arcs_survive_booleans() {
    // A pad with four rounded corners, tagged per element, minus a slot: the remaining
    // corner arcs come back as single arcs.
    let pad = rounded_rect(Point::new(0, 0), 4 * MM, 2 * MM, 500_000);
    let tol = ArcTol::new(100, Side::Nearest);
    let tagged = pad.to_tagged(tol, &|_, j| j as u64 + 1).unwrap();
    // Elements 1, 3, 5, 7 are the corner arcs (centres listed for the lookup).
    let centres = [
        (2, Point::new(1_500_000, -500_000)),
        (4, Point::new(1_500_000, 500_000)),
        (6, Point::new(-1_500_000, 500_000)),
        (8, Point::new(-1_500_000, -500_000)),
    ];
    let arc_of = |t: u64| centres.iter().find(|c| c.0 == t).map(|c| c.1);
    let slot = Ring::from([
        (-100_000, -2 * MM),
        (100_000, -2 * MM),
        (100_000, 2 * MM),
        (-100_000, 2 * MM),
    ]);
    let res = Boolean::new()
        .subject(&tagged, FillRule::NonZero)
        .clip(&slot, FillRule::NonZero)
        .op(Op::Difference)
        .execute_tagged()
        .unwrap();
    assert_eq!(res.len(), 2);
    for p in &res {
        let contour = arcs_from_tags(&p.outer, &arc_of);
        let arcs = contour
            .iter()
            .filter(|c| matches!(c, Curve::CenterArc { .. }))
            .count();
        assert_eq!(arcs, 2, "{contour:?}");
        // Straight parts stay lines; total element count is small.
        assert!(contour.len() <= 8, "{contour:?}");
    }
}
