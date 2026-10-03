use polyclip::*;

fn sq(x0: i64, y0: i64, x1: i64, y1: i64) -> Ring {
    Ring::from([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])
}

fn area2(ps: &PolygonSet) -> i128 {
    ps.iter().map(|p| p.signed_area2()).sum()
}

fn r(v: &[(i64, i64)]) -> Ring {
    Ring::from(v)
}

#[test]
fn union_two_squares() {
    let u = boolean(
        Op::Union,
        &sq(0, 0, 10, 10),
        &sq(5, 5, 15, 15),
        FillRule::NonZero,
    )
    .unwrap();
    assert_eq!(u.len(), 1);
    assert_eq!(
        u[0].outer,
        r(&[
            (0, 0),
            (10, 0),
            (10, 5),
            (15, 5),
            (15, 15),
            (5, 15),
            (5, 10),
            (0, 10)
        ])
    );
    assert_eq!(area2(&u), 2 * 175);
}

#[test]
fn intersection_difference_xor() {
    let a = sq(0, 0, 10, 10);
    let b = sq(5, 5, 15, 15);
    let i = boolean(Op::Intersection, &a, &b, FillRule::NonZero).unwrap();
    assert_eq!(i, vec![Polygon::from(sq(5, 5, 10, 10))]);
    let d = boolean(Op::Difference, &a, &b, FillRule::NonZero).unwrap();
    assert_eq!(
        d[0].outer,
        r(&[(0, 0), (10, 0), (10, 5), (5, 5), (5, 10), (0, 10)])
    );
    let x = boolean(Op::Xor, &a, &b, FillRule::NonZero).unwrap();
    // Two L shapes touching at (5,10)/(10,5)? They touch at the corners (10,5)... no: the
    // xor is two L-shapes meeting at points (10,5)... check areas only.
    assert_eq!(area2(&x), 2 * 150);
    assert_eq!(x.len(), 2);
}

#[test]
fn hole() {
    let d = boolean(
        Op::Difference,
        &sq(0, 0, 10, 10),
        &sq(2, 2, 8, 8),
        FillRule::NonZero,
    )
    .unwrap();
    assert_eq!(d.len(), 1);
    assert_eq!(d[0].outer, sq(0, 0, 10, 10));
    // Hole clockwise, starting at its minimum vertex.
    assert_eq!(d[0].holes, vec![r(&[(2, 2), (2, 8), (8, 8), (8, 2)])]);
    assert_eq!(d[0].signed_area2(), 2 * (100 - 36));
}

#[test]
fn island_tree() {
    let subj = vec![sq(0, 0, 10, 10), sq(4, 4, 6, 6)];
    let t = Boolean::new()
        .subject(&subj, FillRule::NonZero)
        .clip(&sq(2, 2, 8, 8), FillRule::NonZero)
        .op(Op::Xor)
        .execute_tree()
        .unwrap();
    // outer(0..10) > hole(2..8) > island(4..6) minus ... xor: square 4..6 is in both -> out.
    // subject = big square ∪ small (NonZero: small inside big has winding 2 -> inside).
    // xor with 2..8: region 2..8 removed entirely.
    assert_eq!(t.nodes.len(), 2);
    let subj = vec![sq(0, 0, 10, 10)];
    let t = Boolean::new()
        .subject(&subj, FillRule::NonZero)
        .clip(&vec![sq(2, 2, 8, 8), sq(4, 4, 6, 6)], FillRule::EvenOdd)
        .op(Op::Difference)
        .execute_tree()
        .unwrap();
    assert_eq!(t.nodes.len(), 3);
    assert_eq!(t.roots, vec![0]);
    assert!(!t.nodes[0].is_hole);
    assert!(t.nodes[1].is_hole && t.nodes[1].parent == Some(0));
    assert!(!t.nodes[2].is_hole && t.nodes[2].parent == Some(1));
    let ps = t.to_polygon_set();
    assert_eq!(ps.len(), 2);
    assert_eq!(ps[1].outer, sq(4, 4, 6, 6));
}

#[test]
fn touching_corners() {
    let u = union_all(
        &vec![sq(0, 0, 10, 10), sq(10, 10, 20, 20)],
        FillRule::NonZero,
    )
    .unwrap();
    assert_eq!(u.len(), 2);
    assert_eq!(u[0].outer, sq(0, 0, 10, 10));
    assert_eq!(u[1].outer, sq(10, 10, 20, 20));
}

#[test]
fn hole_touching_outer() {
    // Triangle hole touching the bottom edge at (5,0).
    let d = boolean(
        Op::Difference,
        &sq(0, 0, 10, 10),
        &r(&[(5, 0), (7, 5), (3, 5)]),
        FillRule::NonZero,
    )
    .unwrap();
    assert_eq!(d.len(), 1);
    assert_eq!(d[0].holes.len(), 1);
    assert_eq!(d[0].outer, r(&[(0, 0), (5, 0), (10, 0), (10, 10), (0, 10)]));
    assert_eq!(d[0].holes[0], r(&[(3, 5), (7, 5), (5, 0)]));
}

#[test]
fn shared_edge_merges() {
    let u = union_all(
        &vec![sq(0, 0, 10, 10), sq(10, 0, 20, 10)],
        FillRule::NonZero,
    )
    .unwrap();
    assert_eq!(u, vec![Polygon::from(sq(0, 0, 20, 10))]);
}

#[test]
fn self_intersecting_fill_rules() {
    // Pentagram.
    let star = r(&[(0, 10), (6, -8), (-10, 3), (10, 3), (-6, -8)]);
    let eo = union_all(&star, FillRule::EvenOdd).unwrap();
    let nz = union_all(&star, FillRule::NonZero).unwrap();
    assert_eq!(eo.len(), 5);
    assert_eq!(nz.len(), 1);
    assert!(area2(&nz) > area2(&eo));
    // Orientation of the star is clockwise-ish? Positive/Negative pick one side.
    let pos = union_all(&star, FillRule::Positive).unwrap();
    let neg = union_all(&star, FillRule::Negative).unwrap();
    assert!(pos.is_empty() != neg.is_empty());
}

#[test]
fn degenerate_inputs() {
    let empty: Vec<Ring> = vec![];
    assert!(union_all(&empty, FillRule::NonZero).unwrap().is_empty());
    assert!(
        union_all(&r(&[(0, 0), (5, 5)]), FillRule::NonZero)
            .unwrap()
            .is_empty()
    );
    assert!(
        union_all(&r(&[(0, 0), (5, 5), (10, 10)]), FillRule::NonZero)
            .unwrap()
            .is_empty()
    );
    assert!(
        union_all(&r(&[(1, 1), (1, 1), (1, 1)]), FillRule::NonZero)
            .unwrap()
            .is_empty()
    );
    // Same square twice in opposite orientations cancels.
    let mut b = sq(0, 0, 4, 4);
    b.reverse_orientation();
    assert!(
        union_all(&vec![sq(0, 0, 4, 4), b], FillRule::NonZero)
            .unwrap()
            .is_empty()
    );
    // Duplicate square.
    assert_eq!(
        union_all(&vec![sq(0, 0, 4, 4), sq(0, 0, 4, 4)], FillRule::NonZero)
            .unwrap()
            .len(),
        1
    );
    assert!(
        union_all(&vec![sq(0, 0, 4, 4), sq(0, 0, 4, 4)], FillRule::EvenOdd)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn out_of_range() {
    let big = sq(0, 0, MAX_COORD + 1, 10);
    assert!(matches!(
        union_all(&big, FillRule::NonZero),
        Err(Error::CoordinateOutOfRange(_))
    ));
    let ok = sq(-MAX_COORD, -MAX_COORD, MAX_COORD, MAX_COORD);
    let x = r(&[
        (-MAX_COORD, -MAX_COORD),
        (MAX_COORD, MAX_COORD),
        (MAX_COORD, -MAX_COORD),
        (-MAX_COORD, MAX_COORD),
    ]);
    let i = boolean(Op::Intersection, &ok, &x, FillRule::EvenOdd).unwrap();
    assert_eq!(i.len(), 2);
}

#[test]
fn tags_preserved() {
    let a = TaggedRing::uniform(sq(0, 0, 10, 10), 1);
    let b = TaggedRing::uniform(sq(5, 5, 15, 15), 2);
    let res = Boolean::new()
        .subject(&a, FillRule::NonZero)
        .clip(&b, FillRule::NonZero)
        .op(Op::Union)
        .execute_tagged()
        .unwrap();
    let o = &res[0].outer;
    // (0,0)->(10,0)->(10,5) from a; (10,5)->(15,5)->(15,15)->(5,15)->(5,10) from b.
    assert_eq!(o.tags, vec![1, 1, 2, 2, 2, 2, 1, 1]);
}

#[test]
fn clip_open_paths() {
    let square = sq(0, 0, 10, 10);
    let path = Path::from([(-5, 5), (5, 5), (5, 15), (5, 20)]);
    let r = clip_paths(&path, &square, FillRule::NonZero).unwrap();
    assert_eq!(
        r.inside_paths(),
        vec![Path::from([(0, 5), (5, 5), (5, 10)])]
    );
    assert_eq!(
        r.outside_paths(),
        vec![
            Path::from([(-5, 5), (0, 5)]),
            Path::from([(5, 10), (5, 15), (5, 20)])
        ]
    );
    // Along the boundary counts as inside.
    let path = Path::from([(-5, 0), (15, 0)]);
    let r = clip_paths(&path, &square, FillRule::NonZero).unwrap();
    assert_eq!(r.inside_paths(), vec![Path::from([(0, 0), (10, 0)])]);
}

#[test]
fn engine_reuse() {
    let mut e = Boolean::new();
    e.add_subject(&sq(0, 0, 10, 10), FillRule::NonZero);
    let a = e.execute().unwrap();
    e.clear();
    e.add_subject(&sq(0, 0, 10, 10), FillRule::NonZero);
    assert_eq!(e.execute().unwrap(), a);
}

#[test]
fn multi_component_containment() {
    // Regression: the second polygon of a set lies inside `a`; the first does not.
    let a = Ring::from([(-969, 876), (-728, -145), (-323, -400), (0, 284)]);
    let set = vec![
        Polygon::from(Ring::from([(-1189, 0), (-1188, -1), (-1188, 0)])),
        Polygon::from(Ring::from([(-950, 864), (-896, 635), (-391, 0)])),
    ];
    assert!(intersects(&a, &set));
    assert!(intersects(&set, &a));
    assert!(distance(&a, &set).unwrap().sq.is_zero());
    assert!(distance_less_than(&set, &a, 1));
}

#[test]
fn valid_input_unchanged() {
    // Regression: a segment passing exactly through the corner of another vertex's pixel
    // must not be snapped when nothing crosses (union of a valid polygon is the identity).
    let r = Ring::from([(-12, 1), (-11, 0), (-11, 1), (0, 1), (0, 8), (-12, 2)]);
    assert_eq!(
        union_all(&r, FillRule::NonZero).unwrap(),
        vec![Polygon::from(r)]
    );
}
