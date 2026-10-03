use polyclip::*;
use proptest::prelude::*;

fn area(ps: &PolygonSet) -> i128 {
    ps.iter().map(|p| p.signed_area2()).sum()
}

fn rect(x0: i64, y0: i64, w: i64, h: i64, ccw: bool) -> Ring {
    let mut r = Ring::from([(x0, y0), (x0 + w, y0), (x0 + w, y0 + h), (x0, y0 + h)]);
    if !ccw {
        r.reverse_orientation();
    }
    r
}

/// Random rectangles: all intersections are integral, so booleans are exact.
fn rects(range: i64, max: usize) -> impl Strategy<Value = Vec<Ring>> {
    prop::collection::vec(
        (
            0..range,
            0..range,
            1..range / 2 + 2,
            1..range / 2 + 2,
            any::<bool>(),
        )
            .prop_map(|(x, y, w, h, o)| rect(x, y, w, h, o)),
        0..max,
    )
}

/// Arbitrary (usually self-intersecting) rings.
fn rings(range: i64, max_rings: usize, max_pts: usize) -> impl Strategy<Value = Vec<Ring>> {
    prop::collection::vec(
        prop::collection::vec((-range..=range, -range..=range), 0..max_pts)
            .prop_map(|v| v.into_iter().map(Point::from).collect::<Ring>()),
        0..max_rings,
    )
}

fn rule() -> impl Strategy<Value = FillRule> {
    prop_oneof![
        Just(FillRule::EvenOdd),
        Just(FillRule::NonZero),
        Just(FillRule::Positive),
        Just(FillRule::Negative)
    ]
}

fn op() -> impl Strategy<Value = Op> {
    prop_oneof![
        Just(Op::Union),
        Just(Op::Intersection),
        Just(Op::Difference),
        Just(Op::Xor)
    ]
}

fn run(op: Op, a: &[Ring], b: &[Ring], ra: FillRule, rb: FillRule) -> PolygonSet {
    Boolean::new()
        .subject(a, ra)
        .clip(b, rb)
        .op(op)
        .execute()
        .unwrap()
}

proptest! {
    #[test]
    fn rect_identities(a in rects(40, 6), b in rects(40, 6), ra in rule(), rb in rule()) {
        let u = run(Op::Union, &a, &b, ra, rb);
        let i = run(Op::Intersection, &a, &b, ra, rb);
        let d = run(Op::Difference, &a, &b, ra, rb);
        let x = run(Op::Xor, &a, &b, ra, rb);
        let na = union_all(&a, ra).unwrap();
        let nb = union_all(&b, rb).unwrap();
        for s in [&u, &i, &d, &x, &na, &nb] {
            prop_assert_eq!(check_canonical(s, true), Ok(()));
        }
        prop_assert_eq!(area(&u) + area(&i), area(&na) + area(&nb));
        prop_assert_eq!(area(&d) + area(&i), area(&na));
        prop_assert_eq!(area(&x), area(&u) - area(&i));
        // (A - B) ∪ (A ∩ B) = A exactly.
        let back = boolean(Op::Union, &d, &i, FillRule::NonZero).unwrap();
        prop_assert_eq!(&back, &na);
        // A ∩ B ⊆ A and ⊆ B; A ⊆ A ∪ B.
        prop_assert!(contains(&na, &i));
        prop_assert!(contains(&nb, &i));
        prop_assert!(contains(&u, &na));
    }

    #[test]
    fn general_valid(a in rings(12, 4, 8), b in rings(12, 4, 8), ra in rule(), rb in rule(), o in op()) {
        let r = run(o, &a, &b, ra, rb);
        prop_assert_eq!(check_canonical(&r, true), Ok(()));
        let t = Boolean::new().subject(&a, ra).clip(&b, rb).op(o).execute_tree().unwrap();
        prop_assert_eq!(t.to_polygon_set(), r);
    }

    #[test]
    fn general_valid_large(a in rings(1 << 40, 3, 7), b in rings(1 << 40, 3, 7), ra in rule(), o in op()) {
        let r = run(o, &a, &b, ra, ra);
        prop_assert_eq!(check_canonical(&r, true), Ok(()));
    }

    #[test]
    fn general_identities_approx(a in rings(1000, 3, 7), b in rings(1000, 3, 7), ra in rule(), rb in rule()) {
        let u = run(Op::Union, &a, &b, ra, rb);
        let i = run(Op::Intersection, &a, &b, ra, rb);
        let na = union_all(&a, ra).unwrap();
        let nb = union_all(&b, rb).unwrap();
        let lhs = area(&u) + area(&i);
        let rhs = area(&na) + area(&nb);
        // Snap rounding moves edges by < 1 unit: bound by total perimeter.
        let per: f64 = a.iter().chain(b.iter()).map(|r| {
            r.edges().map(|(p, q)| (((q.x - p.x) as f64).powi(2) + ((q.y - p.y) as f64).powi(2)).sqrt()).sum::<f64>()
        }).sum();
        prop_assert!(((lhs - rhs) as f64).abs() <= 8.0 * per + 8.0, "{} vs {} per {}", lhs, rhs, per);
    }

    #[test]
    fn permutation_invariant(a in rings(20, 5, 7), rot in 0usize..7, rule in rule()) {
        let r1 = union_all(&a, rule).unwrap();
        let mut b: Vec<Ring> = a.iter().rev().cloned().collect();
        for r in b.iter_mut() {
            if !r.is_empty() {
                let k = rot % r.len();
                r.rotate_left(k);
            }
        }
        let r2 = union_all(&b, rule).unwrap();
        prop_assert_eq!(r1, r2);
    }

    #[test]
    fn open_path_split(
        paths in prop::collection::vec((-20i64..20, -20i64..20, prop::collection::vec((any::<bool>(), -15i64..15), 0..6)), 0..4),
        clip in rects(20, 4),
    ) {
        // Axis-parallel paths against rectangles: no rounding, exact classification.
        let paths: Vec<Path> = paths.into_iter().map(|(x, y, steps)| {
            let mut p = vec![Point::new(x, y)];
            for (horiz, d) in steps {
                let l = *p.last().unwrap();
                p.push(if horiz { Point::new(l.x + d, l.y) } else { Point::new(l.x, l.y + d) });
            }
            p.into_iter().collect()
        }).collect();
        let r = clip_paths(&paths, &clip, FillRule::NonZero).unwrap();
        let region = union_all(&clip, FillRule::NonZero).unwrap();
        let dbl: PolygonSet = region.iter().map(|poly| Polygon {
            outer: poly.outer.iter().map(|q| Point::new(2 * q.x, 2 * q.y)).collect(),
            holes: poly.holes.iter().map(|h| h.iter().map(|q| Point::new(2 * q.x, 2 * q.y)).collect()).collect(),
        }).collect();
        let total_in: i64 = paths.iter().map(|p| p.windows(2).map(|w| (w[1].x - w[0].x).abs() + (w[1].y - w[0].y).abs()).sum::<i64>()).sum();
        let mut total_out = 0;
        for (pieces, inside) in [(&r.inside, true), (&r.outside, false)] {
            for p in pieces {
                prop_assert!(p.points.len() >= 2);
                prop_assert_eq!(p.tags.len(), p.points.len() - 1);
                for w in p.points.windows(2) {
                    prop_assert_ne!(w[0], w[1]);
                    total_out += (w[1].x - w[0].x).abs() + (w[1].y - w[0].y).abs();
                    let m = Point::new(w[0].x + w[1].x, w[0].y + w[1].y);
                    let loc = locate(&dbl, m);
                    if inside {
                        prop_assert_ne!(loc, Location::Outside);
                    } else {
                        prop_assert_eq!(loc, Location::Outside);
                    }
                }
            }
        }
        prop_assert_eq!(total_in, total_out);
    }
}
