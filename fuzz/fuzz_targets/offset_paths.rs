//! Open-path offsetting (`offset_paths`, `offset_paths_tree`) with every join and end cap.
//!
//! Checks: negative deltas and out-of-range input are rejected; successful output is
//! canonical; the tree form agrees; with a round or square cap and a delta of at least 2,
//! every path vertex is covered by the result.
#![no_main]

use libfuzzer_sys::fuzz_target;
use polyclip::*;
use polyclip_fuzz::{Gen, all_in_range};

fuzz_target!(|data: &[u8]| {
    let mut g = Gen::new(data);
    let paths = g.paths(4, 10);
    let delta = g.delta();
    let join = g.join();
    let cap = g.cap();
    let tol = g.arc_tol_for(delta.saturating_abs());
    polyclip_fuzz::dump(
        "(paths, delta, join, cap, tol)",
        &(&paths, delta, join, cap, tol),
    );
    let in_range = all_in_range(paths.iter().flat_map(|p| p.iter()));

    let out = match offset_paths(&paths, delta, join, cap, tol) {
        Ok(out) => {
            assert!(in_range, "out-of-range input accepted");
            assert!(delta >= 0, "negative path offset accepted");
            out
        }
        Err(Error::CoordinateOutOfRange(_)) | Err(Error::InvalidParameter(_)) => return,
        Err(e) => {
            assert!(!in_range, "unexpected error {e:?}");
            return;
        }
    };
    if let Err(e) = check_canonical(&out, true) {
        panic!("non-canonical path offset output: {e:?}\n{out:?}");
    }
    let tree = offset_paths_tree(&paths, delta, join, cap, tol).expect("tree form");
    assert_eq!(tree.to_polygon_set(), out, "offset_paths_tree disagrees");
    // Covered vertices: path ends under a round/square cap, and interior vertices under
    // any join except bevel (a bevel at a 180° reversal is a flat end through the
    // vertex, which rounding may leave just outside).
    if delta >= 2 {
        for path in &paths {
            let mut v: Vec<Point> = path.to_vec();
            v.dedup();
            let n = v.len();
            for (i, pt) in v.iter().enumerate() {
                let is_end = cap != EndCap::Joined && (i == 0 || i + 1 == n);
                let covered = if is_end {
                    matches!(cap, EndCap::Round | EndCap::Square)
                } else {
                    !matches!(join, Join::Bevel)
                };
                if covered {
                    assert_ne!(
                        locate(&out, *pt),
                        Location::Outside,
                        "path vertex {pt:?} not covered"
                    );
                }
            }
        }
    }
});
