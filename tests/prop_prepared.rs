//! Queries on a `Prepared` geometry must give exactly what the free functions give (closest
//! points included), for valid and invalid geometries of every kind.

use polyclip::*;
use proptest::prelude::*;

#[path = "common/corpus.rs"]
mod corpus;

fn pts(range: i64, n: std::ops::Range<usize>) -> impl Strategy<Value = Vec<Point>> {
    prop::collection::vec((-range..=range, -range..=range), n)
        .prop_map(|v| v.into_iter().map(Point::from).collect())
}

/// Something to query against.
#[derive(Clone, Debug)]
enum Query {
    Ring(Ring),
    Polygon(Polygon),
    Set(Vec<Polygon>),
    Path(Path),
    Point(Point),
    Segment(Segment),
}

fn query(range: i64) -> impl Strategy<Value = Query> {
    prop_oneof![
        3 => pts(range, 3..7).prop_map(|v| Query::Ring(Ring(v))),
        2 => (-range..range, -range..range, 0..range / 4 + 2, 0..range / 4 + 2)
            .prop_map(|(x, y, w, h)| Query::Ring(Ring::from([(x, y), (x + w, y), (x + w, y + h), (x, y + h)]))),
        1 => (pts(range, 3..6), pts(range, 3..5)).prop_map(|(o, h)| Query::Polygon(Polygon::new(Ring(o), vec![Ring(h)]))),
        1 => prop::collection::vec(pts(range, 3..6), 1..4).prop_map(|rs| {
            let rings: Vec<Ring> = rs.into_iter().map(Ring).collect();
            Query::Set(union_all(&rings, FillRule::NonZero).unwrap())
        }),
        2 => pts(range, 1..6).prop_map(|v| Query::Path(Path(v))),
        1 => (-range..=range, -range..=range).prop_map(|(x, y)| Query::Point(Point::new(x, y))),
        1 => pts(range, 2..3).prop_map(|v| Query::Segment(Segment::new(v[0], v[1]))),
    ]
}

fn same_closest(a: Option<Closest>, b: Option<Closest>) -> bool {
    match (a, b) {
        (None, None) => true,
        (Some(a), Some(b)) => {
            a.sq == b.sq
                && a.a.x.to_bits() == b.a.x.to_bits()
                && a.a.y.to_bits() == b.a.y.to_bits()
                && a.b.x.to_bits() == b.b.x.to_bits()
                && a.b.y.to_bits() == b.b.y.to_bits()
        }
        _ => false,
    }
}

fn check_query<G: Preparable + ?Sized, B: Geometry + ?Sized>(
    g: &G,
    p: &Prepared<'_, G>,
    b: &B,
    ds: &[i64],
) -> std::result::Result<(), TestCaseError> {
    prop_assert_eq!(p.intersects(b), intersects(g, b), "intersects");
    prop_assert_eq!(p.contains(b), contains(g, b), "contains");
    for &d in ds {
        prop_assert_eq!(
            p.distance_less_than(b, d),
            distance_less_than(g, b, d),
            "distance_less_than {}",
            d
        );
    }
    let (x, y) = (p.distance(b), distance(g, b));
    prop_assert!(same_closest(x, y), "distance {:?} vs {:?}", x, y);
    Ok(())
}

fn check<G: Preparable + ?Sized>(
    g: &G,
    qs: &[Query],
    probes: &[Point],
    ds: &[i64],
) -> std::result::Result<(), TestCaseError> {
    let p = Prepared::new(g);
    for &q in probes {
        prop_assert_eq!(p.locate(q), locate(g, q), "locate {:?}", q);
    }
    for q in qs {
        match q {
            Query::Ring(b) => check_query(g, &p, b, ds)?,
            Query::Polygon(b) => check_query(g, &p, b, ds)?,
            Query::Set(b) => check_query(g, &p, b, ds)?,
            Query::Path(b) => check_query(g, &p, b, ds)?,
            Query::Point(b) => check_query(g, &p, b, ds)?,
            Query::Segment(b) => check_query(g, &p, b, ds)?,
        }
        // Vertices of the query as probes too (on-boundary cases).
        let mut v = Vec::new();
        match q {
            Query::Ring(b) => v.extend(b.0.iter().copied()),
            Query::Path(b) => v.extend(b.0.iter().copied()),
            Query::Point(b) => v.push(*b),
            _ => {}
        }
        for q in v {
            prop_assert_eq!(p.locate(q), locate(g, q));
        }
    }
    // Every vertex of the prepared geometry is on its boundary (or inside, for overlaps).
    let mut vs = Vec::new();
    g.visit_segments(&mut |a, _| vs.push(a));
    for &q in vs.iter().step_by(7) {
        prop_assert_eq!(p.locate(q), locate(g, q));
        let s = Segment::new(q, Point::new(q.x + 1, q.y + 2));
        check_query(g, &p, &s, ds)?;
    }
    Ok(())
}

/// A zone: a frame minus many small rectangles (valid, canonical, many segments).
fn zone(range: i64) -> impl Strategy<Value = Vec<Polygon>> {
    prop::collection::vec(
        (
            -range..range,
            -range..range,
            1..range / 6 + 2,
            1..range / 6 + 2,
        ),
        20..160,
    )
    .prop_map(move |rs| {
        let frame = Ring::from([
            (-range, -range),
            (range, -range),
            (range, range),
            (-range, range),
        ]);
        let holes: Vec<Ring> = rs
            .iter()
            .map(|&(x, y, w, h)| Ring::from([(x, y), (x + w, y), (x + w, y + h), (x, y + h)]))
            .collect();
        boolean(Op::Difference, &frame, &holes, FillRule::NonZero).unwrap()
    })
}

fn probes(range: i64) -> impl Strategy<Value = Vec<Point>> {
    pts(range + 2, 10..40)
}

/// 64 cases (16 in debug builds) unless `PROPTEST_CASES` says otherwise.
fn cases() -> u32 {
    std::env::var("PROPTEST_CASES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(if cfg!(debug_assertions) { 16 } else { 64 })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(cases()))]

    #[test]
    fn prepared_zone(z in zone(60), qs in prop::collection::vec(query(70), 1..8), pr in probes(60), ds in prop::collection::vec(-1i64..20, 1..4)) {
        check(&z, &qs, &pr, &ds)?;
        for p in z.iter().take(3) {
            check(p, &qs, &pr, &ds)?;
            check(&p.outer, &qs, &pr, &ds)?;
        }
        let tree = Boolean::new().subject(&z, FillRule::NonZero).execute_tree().unwrap();
        check(&tree, &qs, &pr, &ds)?;
    }

    #[test]
    fn prepared_raw(rs in prop::collection::vec(pts(40, 3..200), 1..6), qs in prop::collection::vec(query(45), 1..8), pr in probes(40), ds in prop::collection::vec(-1i64..12, 1..4)) {
        // Self-intersecting, overlapping, any orientation.
        let rings: Vec<Ring> = rs.iter().cloned().map(Ring).collect();
        check(&rings[0], &qs, &pr, &ds)?;
        let poly = Polygon::new(rings[0].clone(), rings[1..].to_vec());
        check(&poly, &qs, &pr, &ds)?;
        let set: Vec<Polygon> = rings.iter().cloned().map(Polygon::from).collect();
        check(&set, &qs, &pr, &ds)?;
        let tree = PolyTree {
            nodes: rings
                .iter()
                .enumerate()
                .map(|(i, r)| PolyNode { ring: r.clone(), tags: vec![0; r.len()], is_hole: i % 2 == 1, parent: None, children: Vec::new() })
                .collect(),
            roots: Vec::new(),
        };
        check(&tree, &qs, &pr, &ds)?;
        let path = Path(rs.concat());
        check(&path, &qs, &pr, &ds)?;
    }

    #[test]
    fn prepared_degenerate(rs in prop::collection::vec(pts(3, 0..40), 1..4), qs in prop::collection::vec(query(5), 1..10), pr in probes(4), ds in prop::collection::vec(-1i64..4, 1..3)) {
        let rings: Vec<Ring> = rs.iter().cloned().map(Ring).collect();
        check(&rings[0], &qs, &pr, &ds)?;
        let poly = Polygon::new(rings[0].clone(), rings[1..].to_vec());
        check(&poly, &qs, &pr, &ds)?;
        let set: Vec<Polygon> = rings.iter().cloned().map(Polygon::from).collect();
        check(&set, &qs, &pr, &ds)?;
        check(&Path(rs.concat()), &qs, &pr, &ds)?;
    }
}

/// A regular grid of holes (many equal distances), at least 512 segments.
fn grid_zone(sizes: &[(i64, i64)]) -> Vec<Polygon> {
    let frame = Ring::from([(0, 0), (480, 0), (480, 480), (0, 480)]);
    let holes: Vec<Ring> = sizes
        .iter()
        .enumerate()
        .map(|(k, &(w, h))| {
            let (x, y) = (8 + 40 * (k as i64 % 12), 8 + 40 * (k as i64 / 12));
            Ring::from([(x, y), (x + w, y), (x + w, y + h), (x, y + h)])
        })
        .collect();
    boolean(Op::Difference, &frame, &holes, FillRule::NonZero).unwrap()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(cases()))]

    #[test]
    fn prepared_distance_ties(sizes in prop::collection::vec((4i64..30, 4i64..30), 144), qs in prop::collection::vec(query(12), 1..12), offs in prop::collection::vec((-30i64..500, -30i64..500), 12), ds in prop::collection::vec(1i64..20, 1..3)) {
        let z = grid_zone(&sizes);
        let p = Prepared::new(&z);
        for (q, &(dx, dy)) in qs.iter().zip(&offs) {
            let mv = |v: &[Point]| -> Vec<Point> { v.iter().map(|p| Point::new(p.x + dx, p.y + dy)).collect() };
            match q {
                Query::Ring(b) => check_query(&z, &p, &Ring(mv(&b.0)), &ds)?,
                Query::Path(b) => check_query(&z, &p, &Path(mv(&b.0)), &ds)?,
                Query::Point(b) => check_query(&z, &p, &Point::new(b.x + dx, b.y + dy), &ds)?,
                Query::Segment(b) => {
                    let v = mv(&[b.a, b.b]);
                    check_query(&z, &p, &Segment::new(v[0], v[1]), &ds)?
                }
                Query::Polygon(b) => {
                    let b = Polygon::new(Ring(mv(&b.outer.0)), b.holes.iter().map(|h| Ring(mv(&h.0))).collect());
                    check_query(&z, &p, &b, &ds)?
                }
                Query::Set(b) => {
                    let b: Vec<Polygon> = b.iter().map(|b| Polygon::new(Ring(mv(&b.outer.0)), b.holes.iter().map(|h| Ring(mv(&h.0))).collect())).collect();
                    check_query(&z, &p, &b, &ds)?
                }
            }
        }
    }
}

#[test]
fn prepared_corpus() {
    let fill = corpus::load("cadlab_gnd_in1.pclp");
    let p = Prepared::new(&fill);
    let b = fill.bbox().unwrap();
    let mut s = 5u64;
    let mut rnd = |m: i64| {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((s >> 33) as i64).rem_euclid(m)
    };
    // Pads around hole vertices (touching, near, overlapping the boundary) and at random.
    let vs: Vec<Point> = fill[0].holes.iter().map(|h| h.0[0]).collect();
    // The free functions are slow on the whole pour (that is the point): fewer in debug.
    for k in 0..if cfg!(debug_assertions) { 25 } else { 150 } {
        let c = if k % 3 == 0 {
            Point::new(
                b.min.x + rnd(b.width().max(1)),
                b.min.y + rnd(b.height().max(1)),
            )
        } else {
            let v = vs[rnd(vs.len() as i64) as usize];
            Point::new(v.x + rnd(400_000) - 200_000, v.y + rnd(400_000) - 200_000)
        };
        let (w, h) = (50_000 + rnd(500_000), 50_000 + rnd(500_000));
        let pad = Ring::from([
            (c.x, c.y),
            (c.x + w, c.y),
            (c.x + w, c.y + h),
            (c.x, c.y + h),
        ]);
        assert_eq!(p.intersects(&pad), intersects(&fill, &pad));
        assert_eq!(
            p.distance_less_than(&pad, 100_000),
            distance_less_than(&fill, &pad, 100_000)
        );
        assert_eq!(p.locate(c), locate(&fill, c));
        assert_eq!(p.contains(&pad), contains(&fill, &pad), "{pad:?}");
        let track = Path::from([(c.x, c.y), (c.x + w, c.y + w), (c.x + w, c.y + w + h)]);
        assert_eq!(p.contains(&track), contains(&fill, &track));
        if k % 5 == 0 {
            assert!(same_closest(p.distance(&pad), distance(&fill, &pad)));
            assert!(same_closest(p.distance(&track), distance(&fill, &track)));
        }
    }
}
