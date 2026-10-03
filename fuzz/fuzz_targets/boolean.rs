//! Boolean operations (`Boolean`, `boolean`, `union_all`) with every operation and fill
//! rule combination.
//!
//! Checks: out-of-range input is rejected with `Err` (and in-range input never is); the
//! output is canonical (`check_canonical`); `execute`, `execute_tree` and
//! `execute_tagged` agree; re-normalizing the output with `union_all` stays canonical
//! and keeps the area up to snap-rounding slack.
#![no_main]

use libfuzzer_sys::fuzz_target;
use polyclip::*;
use polyclip_fuzz::{Gen, all_in_range, perimeter, ring_points};

fuzz_target!(|data: &[u8]| {
    let mut g = Gen::new(data);
    let a = g.rings(5, 12);
    let b = g.rings(5, 12);
    let (ra, rb) = (g.fill_rule(), g.fill_rule());
    let op = g.op();
    let keep = g.int(0, 3) == 0;
    polyclip_fuzz::dump("(a, b, ra, rb, op, keep)", &(&a, &b, ra, rb, op, keep));
    let in_range = all_in_range(ring_points(&a).chain(ring_points(&b)));

    let mut eng = Boolean::new();
    eng.add_subject(&a, ra)
        .add_clip(&b, rb)
        .set_op(op)
        .set_keep_collinear(keep);
    let out = match (eng.execute(), in_range) {
        (Ok(out), true) => out,
        (Err(_), false) => return,
        (Ok(_), false) => panic!("out-of-range input accepted"),
        (Err(e), true) => panic!("in-range input rejected: {e:?}"),
    };
    if let Err(e) = check_canonical(&out, !keep) {
        panic!("non-canonical output: {e:?}\n{out:?}");
    }
    // The engine is reusable: running again gives the same result, also as a tree.
    let tree = eng.execute_tree().expect("second run");
    assert_eq!(tree.to_polygon_set(), out, "execute_tree disagrees");
    let tagged = eng.execute_tagged().expect("tagged run");
    assert_eq!(tagged.len(), out.len());
    for (t, p) in tagged.iter().zip(&out) {
        assert_eq!(t.outer.points, p.outer.0);
        assert_eq!(t.outer.tags.len(), t.outer.points.len(), "one tag per edge");
    }
    // Same through the free function when both rules agree.
    if ra == rb && !keep {
        assert_eq!(boolean(op, &a, &b, ra).expect("boolean()"), out);
    }
    // Re-normalizing canonical output (any rule: its windings are 0/1) stays canonical
    // and keeps the area up to snap-rounding slack. It is *not* always the identity:
    // standard snap rounding is not idempotent (an output edge may pass through the unit
    // pixel of another output vertex and get re-routed there), see
    // tests/fuzz_regressions.rs `canonical_output_not_fixed_point`.
    if !keep {
        let rule = match g.fill_rule() {
            FillRule::Negative => FillRule::Positive,
            r => r,
        };
        let again = union_all(&out, rule).expect("re-normalize");
        if let Err(e) = check_canonical(&again, true) {
            panic!("re-normalized output not canonical: {e:?}");
        }
        let slack = 2.0 * perimeter(&out) + 4.0;
        let diff = (area2(&again) - area2(&out)).unsigned_abs() as f64;
        assert!(diff <= slack, "re-normalizing moved area by {diff}/2");
    }
});
