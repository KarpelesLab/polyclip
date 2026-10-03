//! Compile-time checks that the public types are `Send + Sync` (no global state, engines
//! can be moved between and shared across threads).

use polyclip::*;

fn assert_send_sync<T: Send + Sync>() {}

#[test]
fn public_types_are_send_sync() {
    assert_send_sync::<Point>();
    assert_send_sync::<Ring>();
    assert_send_sync::<Path>();
    assert_send_sync::<Polygon>();
    assert_send_sync::<PolygonSet>();
    assert_send_sync::<PolyTree>();
    assert_send_sync::<TaggedRing>();
    assert_send_sync::<TaggedPath>();
    assert_send_sync::<TaggedPolygon>();
    assert_send_sync::<Boolean>();
    assert_send_sync::<Shape>();
    assert_send_sync::<Circle>();
    assert_send_sync::<ClippedPaths>();
    assert_send_sync::<Closest>();
    assert_send_sync::<SqDist>();
    assert_send_sync::<Triangulation>();
    assert_send_sync::<Trapezoid>();
    assert_send_sync::<Error>();
    assert_send_sync::<ValidityError>();
}

#[test]
fn engines_work_across_threads() {
    let sq = |x: i64| Ring::from([(x, 0), (x + 10, 0), (x + 10, 10), (x, 10)]);
    let results: Vec<PolygonSet> = std::thread::scope(|s| {
        let hs: Vec<_> = (0..4)
            .map(|i| {
                s.spawn(move || union_all(&vec![sq(0), sq(5 + i)], FillRule::NonZero).unwrap())
            })
            .collect();
        hs.into_iter().map(|h| h.join().unwrap()).collect()
    });
    for (i, r) in results.iter().enumerate() {
        assert_eq!(r[0].outer.signed_area2(), 2 * 10 * (15 + i as i128));
    }
}
