//! Property tests for triangulation: random valid polygons (normalized by the boolean
//! engine, so holes often touch each other and the outer ring at vertices) must give
//! exact, non-overlapping, counter-clockwise triangulations that contain every boundary
//! edge; arbitrary input must never panic.

use polyclip::predicates::orient;
use polyclip::*;
use proptest::prelude::*;
use std::collections::HashMap;

/// Arbitrary (usually self-intersecting) rings.
fn rings(range: i64, max_rings: usize, max_pts: usize) -> impl Strategy<Value = Vec<Ring>> {
    prop::collection::vec(
        prop::collection::vec((-range..=range, -range..=range), 0..max_pts)
            .prop_map(|v| v.into_iter().map(Point::from).collect::<Ring>()),
        0..max_rings,
    )
}

/// Axis-aligned rectangles on a coarse grid (many shared corners).
fn rects(range: i64, max: usize) -> impl Strategy<Value = Vec<Ring>> {
    prop::collection::vec(
        (0..range, 0..range, 1..range / 2 + 2, 1..range / 2 + 2)
            .prop_map(|(x, y, w, h)| Ring::from([(x, y), (x + w, y), (x + w, y + h), (x, y + h)])),
        0..max,
    )
}

fn rule() -> impl Strategy<Value = FillRule> {
    prop_oneof![
        Just(FillRule::EvenOdd),
        Just(FillRule::NonZero),
        Just(FillRule::Positive),
    ]
}

fn normalize(rings: &[Ring], rule: FillRule, keep_collinear: bool) -> PolygonSet {
    Boolean::new()
        .subject(rings, rule)
        .keep_collinear(keep_collinear)
        .execute()
        .unwrap()
}

/// Canonically oriented directed boundary edges.
fn boundary(polys: &[Polygon]) -> Vec<(Point, Point)> {
    let mut out = Vec::new();
    for poly in polys {
        for (ri, r) in poly.rings().enumerate() {
            let mut v = r.0.clone();
            if (ri == 0) != (r.signed_area2() > 0) {
                v.reverse();
            }
            for i in 0..v.len() {
                out.push((v[i], v[(i + 1) % v.len()]));
            }
        }
    }
    out
}

fn incircle_i128(a: Point, b: Point, c: Point, d: Point) -> i128 {
    let f = |u: Point| ((u.x - d.x) as i128, (u.y - d.y) as i128);
    let ((ax, ay), (bx, by), (cx, cy)) = (f(a), f(b), f(c));
    (ax * ax + ay * ay) * (bx * cy - by * cx)
        + (bx * bx + by * by) * (cx * ay - cy * ax)
        + (cx * cx + cy * cy) * (ax * by - ay * bx)
}

/// Checks all invariants of a triangulation of the valid polygon set `polys`.
fn check(
    polys: &[Polygon],
    t: &Triangulation,
    delaunay: bool,
) -> core::result::Result<(), TestCaseError> {
    prop_assert!(t.vertices.windows(2).all(|w| w[0] < w[1]));
    prop_assert!(t.triangles.windows(2).all(|w| w[0] < w[1]));
    let mut half: HashMap<(u32, u32), u32> = HashMap::new();
    let mut sum: i128 = 0;
    for (i, tri) in t.triangles.iter().enumerate() {
        prop_assert!(tri[0] < tri[1] && tri[0] < tri[2]);
        let [a, b, c] = t.triangle(i);
        let o = orient(a, b, c);
        prop_assert!(o > 0, "degenerate or clockwise triangle {:?}", (a, b, c));
        sum += o;
        for k in 0..3 {
            *half.entry((tri[k], tri[(k + 1) % 3])).or_default() += 1;
        }
    }
    // Canonical polygons: outer counter-clockwise, holes clockwise.
    let area: i128 = polys.iter().map(|p| p.signed_area2()).sum();
    prop_assert_eq!(sum, area, "area");
    prop_assert_eq!(t.area2(), area);
    let mut verts: Vec<Point> = Vec::new();
    let mut bnd: HashMap<(u32, u32), u32> = HashMap::new();
    for (a, b) in boundary(polys) {
        verts.push(a);
        let ia = t.vertices.binary_search(&a);
        let ib = t.vertices.binary_search(&b);
        prop_assert!(ia.is_ok() && ib.is_ok(), "boundary vertex missing");
        *bnd.entry((ia.unwrap() as u32, ib.unwrap() as u32))
            .or_default() += 1;
    }
    verts.sort_unstable();
    verts.dedup();
    prop_assert_eq!(&verts, &t.vertices, "vertex set");
    for (&(a, b), &n) in &bnd {
        prop_assert_eq!(n, 1);
        prop_assert_eq!(
            half.get(&(a, b)),
            Some(&1),
            "boundary edge not a triangle edge"
        );
        prop_assert_eq!(half.get(&(b, a)), None, "triangle outside the polygon");
    }
    // Every interior edge is shared by exactly two triangles in opposite directions.
    // Together with positive orientation and the exact area match, this rules out any
    // overlap (the triangles form a degree-one covering of the polygon).
    for (&(a, b), &n) in &half {
        prop_assert_eq!(n, 1);
        if !bnd.contains_key(&(a, b)) {
            prop_assert_eq!(half.get(&(b, a)), Some(&1), "unmatched interior edge");
        }
    }
    if delaunay {
        let mut opp: HashMap<(u32, u32), u32> = HashMap::new();
        for tri in &t.triangles {
            for k in 0..3 {
                opp.insert((tri[k], tri[(k + 1) % 3]), tri[(k + 2) % 3]);
            }
        }
        let v = |i: u32| t.vertices[i as usize];
        for tri in &t.triangles {
            for k in 0..3 {
                let (a, b, c) = (tri[k], tri[(k + 1) % 3], tri[(k + 2) % 3]);
                if let Some(&d) = opp.get(&(b, a)) {
                    prop_assert!(incircle_i128(v(a), v(b), v(c), v(d)) <= 0, "not Delaunay");
                }
            }
        }
    }
    Ok(())
}

fn check_all(polys: &PolygonSet) -> core::result::Result<(), TestCaseError> {
    for poly in polys {
        let one = core::slice::from_ref(poly);
        let t = triangulate(poly).unwrap();
        check(one, &t, false)?;
        let d = triangulate_delaunay(poly).unwrap();
        // The direct i128 in-circle evaluation used by the check overflows beyond 2^20.
        let small = poly
            .outer
            .iter()
            .all(|p| p.x.abs() < 1 << 20 && p.y.abs() < 1 << 20);
        check(one, &d, small)?;
        prop_assert_eq!(&triangulate(poly).unwrap(), &t, "deterministic");
    }
    let t = triangulate_set(polys).unwrap();
    check(polys, &t, false)?;
    Ok(())
}

proptest! {
    #[test]
    fn random_rings(rs in rings(40, 6, 14), rule in rule(), keep in any::<bool>()) {
        check_all(&normalize(&rs, rule, keep))?;
    }

    #[test]
    fn random_rings_wide(rs in rings(1_000_000, 5, 20), rule in rule()) {
        check_all(&normalize(&rs, rule, false))?;
    }

    #[test]
    fn random_rings_full_range(rs in rings(MAX_COORD, 4, 12), rule in rule()) {
        check_all(&normalize(&rs, rule, false))?;
    }

    #[test]
    fn grid_rects(rs in rects(12, 14), rule in rule(), keep in any::<bool>()) {
        // Even-odd of overlapping rectangles creates holes touching at corners.
        check_all(&normalize(&rs, rule, keep))?;
    }

    #[test]
    fn boxed_difference(rs in rings(30, 8, 8)) {
        // A box minus random shapes: many holes, touching each other and the box.
        let outline = Ring::from([(-30, -30), (30, -30), (30, 30), (-30, 30)]);
        let res = boolean(Op::Difference, &outline, &rs, FillRule::NonZero).unwrap();
        check_all(&res)?;
    }

    #[test]
    fn boxed_rects(rs in rects(16, 16), rule in rule()) {
        // A box minus rectangles: holes touching each other (and the box) at corners.
        let outline = Ring::from([(-1, -1), (30, -1), (30, 30), (-1, 30)]);
        let res = boolean(Op::Difference, &outline, &rs, rule).unwrap();
        check_all(&res)?;
        let res = boolean(Op::Xor, &outline, &rs, rule).unwrap();
        check_all(&res)?;
    }

    #[test]
    fn arbitrary_input_never_panics(
        rs in prop_oneof![rings(20, 5, 10), rings(MAX_COORD, 4, 8)],
        big in any::<bool>(),
    ) {
        let mut rs = rs;
        if big && !rs.is_empty() && !rs[0].is_empty() {
            rs[0][0] = Point::new(MAX_COORD + 1, 0);
        }
        let poly = match rs.split_first() {
            Some((o, h)) => Polygon::new(o.clone(), h.to_vec()),
            None => Polygon::default(),
        };
        for t in [triangulate(&poly), triangulate_delaunay(&poly)].into_iter().flatten() {
            for i in 0..t.triangles.len() {
                let [a, b, c] = t.triangle(i);
                prop_assert!(orient(a, b, c) > 0);
            }
        }
        let _ = triangulate_set(&[poly.clone(), poly]);
    }
}
