//! Offsets (whose final union usually splits into many independent clusters of raw rings)
//! give the same result computed cluster by cluster as in one piece.
//!
//! Uses the per-thread switch to force the one-piece computation.

mod common;

use polyclip::*;

#[test]
fn offsets_clustered_identical() {
    let fill = common::corpus::load("cadlab_gnd_in1.pclp");
    let tol = ArcTol::new(5_000, Side::Inside);
    let run = |mono: bool| {
        set_always_monolithic(mono);
        let opened = opening(&fill, 100_000, tol).unwrap();
        // Slow in debug builds: release only.
        let closed = (!cfg!(debug_assertions)).then(|| closing(&fill, 60_000, tol).unwrap());
        let grown = offset_tree(&fill, 30_000, Join::Miter { limit: 2.0 }, tol).unwrap();
        set_always_monolithic(false);
        (opened, closed, grown)
    };
    let (a, b) = (run(false), run(true));
    assert!(a.0 == b.0, "opening differs");
    assert!(a.1 == b.1, "closing differs");
    assert!(a.2 == b.2, "grown tree differs");
}
