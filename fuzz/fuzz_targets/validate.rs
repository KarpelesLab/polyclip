//! Validity checks (`validate`, `validate_set`, `check_canonical`) on arbitrary polygon
//! sets, including out-of-range coordinates, and on near-valid sets made by perturbing
//! canonical output.
//!
//! Checks: no panic; out-of-range coordinates are reported as invalid;
//! `check_canonical` implies `validate_set`, which implies `validate` of each polygon;
//! `union_all` of a canonical set is canonical; a set
//! accepted by `validate_set` keeps its area under `union_all` up to snap-rounding slack.
#![no_main]

use libfuzzer_sys::fuzz_target;
use polyclip::*;
use polyclip_fuzz::{Gen, all_in_range, perimeter, poly_points};

fn check(set: &[Polygon]) {
    polyclip_fuzz::dump("set", &set);
    let valid = validate_set(set);
    let canon = check_canonical(set, true);
    let canon_c = check_canonical(set, false);
    if !all_in_range(poly_points(set)) {
        assert!(valid.is_err(), "out-of-range set validated");
        return;
    }
    if canon.is_ok() {
        assert!(canon_c.is_ok(), "canonical(true) but not canonical(false)");
    }
    if canon_c.is_ok() {
        assert!(valid.is_ok(), "canonical but not valid");
    }
    if valid.is_ok() {
        for p in set {
            assert_eq!(validate(p), Ok(()), "set valid but member invalid");
        }
        // Holes may have either orientation for `validate_set`: normalize by role.
        let fixed: Vec<Polygon> = set
            .iter()
            .map(|p| {
                let mut p = p.clone();
                if !p.outer.is_ccw() {
                    p.outer.reverse_orientation();
                }
                for h in &mut p.holes {
                    if h.is_ccw() {
                        h.reverse_orientation();
                    }
                }
                p
            })
            .collect();
        // No edge crosses another, but snap rounding still routes an edge through any
        // vertex whose unit pixel it passes, so the area may move by a sliver along the
        // boundary (a thin triangle can even collapse).
        let u = union_all(&fixed, FillRule::NonZero).expect("union of valid set");
        let slack = 2.0 * perimeter(&fixed) + 4.0;
        let diff = (area2(&u) - area2(&fixed)).unsigned_abs() as f64;
        assert!(
            diff <= slack,
            "valid set changed area by {diff}/2 under union_all (slack {slack})"
        );
    }
    if canon.is_ok() {
        // Not necessarily the identity (snap rounding is not idempotent), but canonical.
        let u = union_all(set, FillRule::NonZero).expect("union of canonical set");
        if let Err(e) = check_canonical(&u, true) {
            panic!("union_all of a canonical set is not canonical: {e:?}");
        }
    }
}

fuzz_target!(|data: &[u8]| {
    let mut g = Gen::new(data);
    let raw = g.raw_polygons(3);
    check(&raw);
    for p in &raw {
        let _ = validate(p);
    }
    // Near-valid: canonical output, then a small perturbation.
    let rings = g.rings(4, 10);
    let Ok(mut set) = union_all(&rings, g.fill_rule()) else {
        return;
    };
    check(&set);
    if set.is_empty() {
        return;
    }
    let pi = g.int(0, set.len() as i64 - 1) as usize;
    let nr = 1 + set[pi].holes.len();
    let ri = g.int(0, nr as i64 - 1) as usize;
    let ring = if ri == 0 {
        &mut set[pi].outer
    } else {
        &mut set[pi].holes[ri - 1]
    };
    let vi = g.int(0, ring.len() as i64 - 1) as usize;
    match g.int(0, 5) {
        0 => {
            let d = g.int(-2, 2);
            ring[vi].x = ring[vi].x.saturating_add(d);
        }
        1 => {
            let v = ring[vi];
            ring.insert(vi, v);
        }
        2 => {
            ring.remove(vi);
        }
        3 => ring.reverse_orientation(),
        4 => {
            let k = vi;
            ring.rotate_left(k);
        }
        _ => {
            let p = set[pi].clone();
            set.push(p);
        }
    }
    check(&set);
});
