use polyclip::*;
fn main() {
    let a = vec![Ring::from([
        (-323, -400),
        (-728, -145),
        (-969, 876),
        (0, 284),
    ])];
    let dx = -1188;
    let b = vec![
        Ring::from([(0, 0), (0, -1), (-1, 0)]),
        Ring::from([(292, 635), (238, 864), (797, 0)]),
    ];
    let b: Vec<Ring> = b
        .iter()
        .map(|r| r.iter().map(|p| Point::new(p.x + dx, p.y)).collect())
        .collect();
    let pa = union_all(&a, FillRule::NonZero).unwrap();
    let pb = union_all(&b, FillRule::NonZero).unwrap();
    println!("pa {:?}\npb {:?}", pa, pb);
    let c = distance(&pa, &pb).unwrap();
    println!("dist {} {:?} {:?}", c.sq.distance_f64(), c.a, c.b);
    println!("intersects {}", intersects(&pa, &pb));
    let i = boolean(Op::Intersection, &pa, &pb, FillRule::NonZero).unwrap();
    println!("inter {:?}", i);
}
