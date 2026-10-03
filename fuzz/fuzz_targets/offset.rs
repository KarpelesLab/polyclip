//! Region offsetting (`offset`, `offset_tree`, `opening`, `closing`) with every join type,
//! arc side and tolerance, including degenerate miter limits and huge deltas.
//!
//! Checks: out-of-range input is rejected; in-range input with a valid tolerance and delta
//! either succeeds or fails with `CoordinateOutOfRange` (the offset left the coordinate
//! range) / `InvalidParameter`; successful output is canonical; `offset_tree` agrees with
//! `offset`; growing never loses area and shrinking never gains area relative to the
//! normalized input.
#![no_main]

use libfuzzer_sys::fuzz_target;
use polyclip::*;
use polyclip_fuzz::{Gen, all_in_range, ring_points};

fuzz_target!(|data: &[u8]| {
    let mut g = Gen::new(data);
    let rings = g.rings(4, 10);
    let delta = g.delta();
    let join = g.join();
    let tol = g.arc_tol_for(delta.saturating_abs());
    let mode = g.int(0, 3);
    polyclip_fuzz::dump(
        "(rings, delta, join, tol, mode)",
        &(&rings, delta, join, tol, mode),
    );
    let in_range = all_in_range(ring_points(&rings));

    let res = match mode {
        0 => opening(&rings, delta, tol),
        1 => closing(&rings, delta, tol),
        _ => offset(&rings, delta, join, tol),
    };
    let out = match res {
        Ok(out) => {
            assert!(in_range, "out-of-range input accepted");
            out
        }
        Err(Error::CoordinateOutOfRange(_)) | Err(Error::InvalidParameter(_)) => return,
        Err(e) => {
            assert!(!in_range, "unexpected error {e:?}");
            return;
        }
    };
    if let Err(e) = check_canonical(&out, true) {
        panic!("non-canonical offset output: {e:?}\n{out:?}");
    }
    if mode >= 2 {
        let tree = offset_tree(&rings, delta, join, tol).expect("offset_tree");
        assert_eq!(tree.to_polygon_set(), out, "offset_tree disagrees");
        let norm = union_all(&rings, FillRule::NonZero).expect("normalize");
        let (an, ao) = (area2(&norm), area2(&out));
        if delta > 0 {
            assert!(ao >= an, "growing lost area: {an} -> {ao}");
        } else if delta < 0 {
            assert!(ao <= an, "shrinking gained area: {an} -> {ao}");
        } else {
            assert_eq!(out, norm, "zero offset is not normalization");
        }
    }
});
