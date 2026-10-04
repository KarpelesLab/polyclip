//! Real-world corpus: zone fills dumped from cadlab boards (`testdata/`).

mod common;

use polyclip::*;

/// cadlab's synthetic large board, GND pour on In1.Cu after clearances, before the
/// minimum-width opening: one polygon with 2153 holes, 242 316 vertices (POLYGON_LIB.md §8).
#[test]
fn cadlab_gnd_in1_opening() {
    let fill = common::corpus::load("cadlab_gnd_in1.pclp");
    assert_eq!(fill.len(), 1);
    assert_eq!(fill[0].holes.len(), 2153);
    assert_eq!(
        fill.iter().map(|p| p.vertex_count()).sum::<usize>(),
        242_316
    );
    assert_eq!(check_canonical(&fill, true), Ok(()));
    // Canonical input is a fixed point of normalization.
    assert_eq!(union_all(&fill, FillRule::NonZero).unwrap(), fill);
    let tol = ArcTol::new(5_000, Side::Inside);
    let shrunk = offset(&fill, -100_000, Join::Round, tol).unwrap();
    assert_eq!(check_canonical(&shrunk, true), Ok(()));
    let opened = opening(&fill, 100_000, tol).unwrap();
    assert_eq!(check_canonical(&opened, true), Ok(()));
    assert_eq!(opened, offset(&shrunk, 100_000, Join::Round, tol).unwrap());
    assert!(area2(&opened) <= area2(&fill));
    // Thermal spokes: union with 40 small rectangles.
    let spokes: Vec<Ring> = (0..40)
        .map(|k| {
            let (x, y) = (10_000_000 + k * 3_000_000, 50_000_000);
            Ring::from([
                (x, y),
                (x + 250_000, y),
                (x + 250_000, y + 1_000_000),
                (x, y + 1_000_000),
            ])
        })
        .collect();
    let with_spokes = Boolean::new()
        .subject(&opened, FillRule::NonZero)
        .subject(&spokes, FillRule::NonZero)
        .execute()
        .unwrap();
    assert_eq!(check_canonical(&with_spokes, true), Ok(()));
    assert!(area2(&with_spokes) >= area2(&opened));
}
