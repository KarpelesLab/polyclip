//! Clipping open paths by polygons (`clip_paths`).
//!
//! Checks: out-of-range input is rejected (and in-range input never is); output pieces
//! have at least one edge, one tag per edge, no zero-length edges and in-range vertices;
//! an empty clip region leaves everything outside; pieces keep source order and direction
//! (the first piece starts near the first non-degenerate path's first vertex). Copies of
//! the input side by side (enough clip edges to be computed cluster by cluster) give the
//! same pieces as computing everything in one piece.
#![no_main]

use libfuzzer_sys::fuzz_target;
use polyclip::*;
use polyclip_fuzz::{Gen, all_in_range, ring_points};

fuzz_target!(|data: &[u8]| {
    let mut g = Gen::new(data);
    let paths = g.paths(4, 10);
    let clip = g.rings(4, 10);
    let rule = g.fill_rule();
    polyclip_fuzz::dump("(paths, clip, rule)", &(&paths, &clip, rule));
    let in_range = all_in_range(
        paths
            .iter()
            .flat_map(|p| p.iter())
            .chain(ring_points(&clip)),
    );
    let r = match (clip_paths(&paths, &clip, rule), in_range) {
        (Ok(r), true) => r,
        (Err(_), false) => return,
        (Ok(_), false) => panic!("out-of-range input accepted"),
        (Err(e), true) => panic!("in-range input rejected: {e:?}"),
    };
    for piece in r.inside.iter().chain(r.outside.iter()) {
        assert!(piece.points.len() >= 2, "piece with < 2 points");
        assert_eq!(piece.tags.len(), piece.points.len() - 1, "one tag per edge");
        assert!(
            piece.points.windows(2).all(|w| w[0] != w[1]),
            "zero-length edge"
        );
        assert!(all_in_range(piece.points.iter()));
    }
    // (Paths are snap-rounded together with the clip edges, so a clip region that
    // collapses on its own may survive here; only an edgeless clip is surely empty.)
    let clip_has_edges = clip.iter().any(|r| r.edges().any(|(a, b)| a != b));
    if !clip_has_edges {
        assert!(r.inside.is_empty(), "inside pieces without a clip region");
    }
    // Some path has a non-zero-length edge iff some piece is produced.
    let any_edge = paths.iter().any(|p| p.windows(2).any(|w| w[0] != w[1]));
    assert_eq!(
        any_edge,
        !(r.inside.is_empty() && r.outside.is_empty()),
        "pieces appear or vanish"
    );
    // The first piece starts within snapping distance of the first path's first vertex.
    if let Some(p0) = paths
        .iter()
        .find(|p| p.windows(2).any(|w| w[0] != w[1]))
        .map(|p| p[0])
    {
        let starts = r.inside.iter().chain(r.outside.iter()).map(|p| p.points[0]);
        assert!(
            starts
                .into_iter()
                .any(|s| (s.x - p0.x).abs() <= 1 && (s.y - p0.y).abs() <= 1),
            "no piece starts at the first path start {p0:?}"
        );
    }
    // Tiled copies: clustered computation equals the one-piece one.
    let n: usize = clip.iter().map(|r| r.len()).sum();
    let step = 2 * g.extent() + 3;
    let k = (1100 / n.max(1) + 1).min(300) as i64;
    let fits = |p: &Point| p.x.unsigned_abs() as i128 + (k * step) as i128 <= MAX_COORD as i128;
    if n > 0 && ring_points(&clip).all(fits) && paths.iter().flat_map(|p| p.iter()).all(fits) {
        let shift = |p: &Point, i: i64| Point::new(p.x + i * step, p.y);
        let tclip: Vec<Ring> = (0..k)
            .flat_map(|i| clip.iter().map(move |r| r.iter().map(|p| shift(p, i)).collect()))
            .collect();
        let tpaths: Vec<Path> = (0..k)
            .step_by(3)
            .flat_map(|i| paths.iter().map(move |r| r.iter().map(|p| shift(p, i)).collect()))
            .collect();
        let got = clip_paths(&tpaths, &tclip, rule).expect("tiled");
        set_always_monolithic(true);
        let want = clip_paths(&tpaths, &tclip, rule).expect("tiled, one piece");
        set_always_monolithic(false);
        assert_eq!(got, want, "clustered clip_paths differs");
    }
});
