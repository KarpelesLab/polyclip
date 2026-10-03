//! Library bugs found by the fuzz targets in `fuzz/`, decoded into readable test cases.
//!
//! Each test is `#[ignore]`d while the bug is open (the fuzz targets carry a matching
//! guard so they can keep searching for new bugs); remove the `#[ignore]` and the guard
//! once fixed. Run them with `cargo test --test fuzz_regressions -- --ignored`.

use polyclip::*;

fn p(x: i64, y: i64) -> Point {
    Point::new(x, y)
}

/// Found by `fuzz/boolean`: the union of a unit square and a self-overlapping ring trips
/// `debug_assert!(.., "sweep: non-noded input")` in src/sweep.rs (the snap-rounded
/// arrangement handed to the sweep still has an unresolved contact). Each operand alone
/// is fine. In release builds the assertion is compiled out; see the companion check
/// below for whether the output is still valid there.
#[test]
#[ignore = "open bug: snap-rounded arrangement not fully noded (sweep debug_assert)"]
fn union_square_and_self_overlapping_ring_non_noded() {
    let a = vec![Ring(vec![p(3, 1), p(3, 0), p(2, 0), p(2, 1)])];
    let b = vec![Ring(vec![p(3, 1), p(1, 3), p(1, 1), p(1, 5), p(0, 5)])];
    for keep in [false, true] {
        let out = Boolean::new()
            .subject(&a, FillRule::EvenOdd)
            .clip(&b, FillRule::EvenOdd)
            .op(Op::Union)
            .keep_collinear(keep)
            .execute()
            .unwrap();
        assert_eq!(check_canonical(&out, !keep), Ok(()));
    }
}

/// Same `sweep: non-noded input` debug assertion, reached from the plainest CAD input:
/// morphological closing of a 1x2 rectangle by 1 unit (found by `fuzz/offset`; the
/// `offset_paths` and `arc` targets hit it too).
#[test]
#[ignore = "open bug: snap-rounded arrangement not fully noded (sweep debug_assert)"]
fn closing_small_rectangle_non_noded() {
    let r = vec![Ring(vec![p(2, 1), p(2, 3), p(1, 3), p(1, 1)])];
    let out = closing(&r, 1, ArcTol::new(9, Side::Outside)).unwrap();
    assert_eq!(check_canonical(&out, true), Ok(()));
}

/// Found by `fuzz/offset_paths`: same assertion from a mitered closed stroke.
#[test]
#[ignore = "open bug: snap-rounded arrangement not fully noded (sweep debug_assert)"]
fn joined_path_offset_non_noded() {
    let paths = vec![Path(vec![p(3, 3), p(0, 0), p(1, 3), p(1, 2)])];
    let out = offset_paths(
        &paths,
        1,
        Join::Miter { limit: 0.0 },
        EndCap::Joined,
        ArcTol::new(1, Side::Outside),
    )
    .unwrap();
    assert_eq!(check_canonical(&out, true), Ok(()));
}

/// Found by `fuzz/boolean`. Not a validity bug (both results are canonical), but worth
/// knowing: canonical output is not a fixed point of `union_all`. Normalizing the ring
/// below yields two triangles; the edge `(3,1)-(1,0)` of the second one is a snapped
/// fragment that passes through the unit pixel of the vertex `(2,0)`, so normalizing the
/// output again re-routes it through `(2,0)` and the triangle (area 1/2) vanishes.
/// Standard snap rounding is not idempotent; iterated snap rounding (or a final
/// "no edge through a foreign hot pixel" pass) would make it so.
#[test]
#[ignore = "design question: snap-rounded output is not idempotent under union_all"]
fn canonical_output_not_fixed_point() {
    let a = vec![Ring(vec![
        p(3, 1),
        p(3, 1),
        p(3, 1),
        p(2, 0),
        p(0, 1),
        p(0, 0),
    ])];
    let once = union_all(&a, FillRule::NonZero).unwrap();
    assert_eq!(check_canonical(&once, true), Ok(()));
    assert_eq!(once.len(), 2);
    let twice = union_all(&once, FillRule::NonZero).unwrap();
    assert_eq!(twice, once);
}

/// Found by `fuzz/clip_paths`: a single-vertex path outside `±MAX_COORD` is accepted
/// (the range check only looks at vertices of non-degenerate edges), where every other
/// operation returns `Error::CoordinateOutOfRange`.
#[test]
#[ignore = "open bug: clip_paths accepts an out-of-range single-vertex path"]
fn clip_paths_single_vertex_out_of_range() {
    let paths = vec![Path(vec![p(MAX_COORD + 1, MAX_COORD + 1)])];
    let clip: Vec<Ring> = vec![];
    assert!(matches!(
        clip_paths(&paths, &clip, FillRule::EvenOdd),
        Err(Error::CoordinateOutOfRange(_))
    ));
}

/// Found by `fuzz/offset`: `opening` / `closing` compute `-d.abs()` before validating
/// `d`, which overflows (panics in debug builds) for `d == i64::MIN`.
#[test]
#[ignore = "open bug: opening/closing overflow on i64::MIN"]
fn opening_closing_i64_min() {
    let r: Vec<Ring> = vec![];
    let tol = ArcTol::new(1, Side::Outside);
    assert!(opening(&r, i64::MIN, tol).is_err());
    assert!(closing(&r, i64::MIN, tol).is_err());
}

/// Found by `fuzz/arc`: `Shape::to_polygon` computes the contour orientation in exact
/// integer arithmetic before range-checking the input, so out-of-range points overflow
/// (`attempt to subtract with overflow` in `predicates::orient`) instead of returning
/// `Error::CoordinateOutOfRange`.
#[test]
#[ignore = "open bug: Shape::to_polygon overflows on out-of-range points"]
fn shape_out_of_range_overflows() {
    let m = MAX_COORD + 1;
    let shape = Shape::new(
        vec![
            Curve::Arc {
                mid: p(3, m),
                end: p(3, -i64::MAX),
            },
            Curve::Line(p(m, m)),
            Curve::Line(p(m, m)),
            Curve::Line(p(m, m)),
        ],
        vec![],
    );
    assert!(shape.to_polygon(ArcTol::new(1, Side::Outside)).is_err());
}
