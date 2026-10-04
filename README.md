# polyclip

[![CI](https://github.com/KarpelesLab/polyclip/actions/workflows/ci.yml/badge.svg)](https://github.com/KarpelesLab/polyclip/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/polyclip.svg)](https://crates.io/crates/polyclip)
[![docs.rs](https://img.shields.io/docsrs/polyclip)](https://docs.rs/polyclip)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Exact integer 2D polygon geometry for Rust: boolean operations, offsetting, arc
approximation, distance queries, fracturing, simplification and more — with exact
predicates, guaranteed-valid output after rounding, and deterministic canonical results.

Built for the cadlab PCB tool (zone fills, DRC, silkscreen clipping, pad shapes, Gerber
output), usable as a standalone crate.

## Features

- **Booleans**: union, intersection, difference, xor; even-odd / non-zero / positive /
  negative fill rules per operand; N-ary union in one sweep; open-path clipping
  (inside/outside pieces); `PolygonSet` or full `PolyTree` output.
- **Offsetting**: polygons with holes (round, miter, bevel, square joins), open paths
  (round, square, butt caps; closed loops), opening/closing for minimum-width enforcement,
  curved shapes.
- **Incremental zone refill**: `ZoneFill` adds/removes/moves obstacles by id and updates
  the fill in well under a millisecond, with results identical to a full recompute.
- **Curved booleans**: `curved_boolean` on `Shape`s keeps arcs as arcs (for Gerber/IPC-2581
  output) with a guaranteed error bound and side.
- **Arcs**: lines, three-point arcs, centre arcs and circles approximated within a
  tolerance on a selectable side (`Outside`, `Inside`, `Nearest`) so clearances are never
  under-estimated.
- **Provenance tags**: every edge can carry a `u64` tag that survives booleans and offsets
  (reconstruct arcs for Gerber/IPC-2581, explain DRC results).
- **Queries**: exact area, centroid, bounding box, point location, `intersects`,
  `contains`, exact minimum distance with closest points, and a fast
  `distance_less_than` for DRC.
- **Utilities**: validity check with reasons, fracture (holes joined by zero-width cut-ins
  for Gerber regions), topology-preserving simplification, convex hull, Minkowski sum.

## Guarantees

- Integer `i64` coordinates, supported range `±2^40` (out-of-range input is an error).
- Exact predicates (`i128` / wide integer arithmetic); floats only for constructing new
  vertices and as filters in front of exact fallbacks.
- Snap rounding: outputs are valid after rounding (simple rings, no crossings, correct
  nesting), vertices move by at most `√2/2`, and input that needs no rounding passes
  through unchanged.
- Canonical output: outer rings CCW, holes CW, rings start at their smallest vertex,
  collinear vertices removed, sorted polygons and holes; bit-identical across platforms.
- Never panics; no `unsafe`; no global state.

## Example

```rust
use polyclip::*;

fn main() -> polyclip::Result<()> {
    let zone = Ring::from([(0, 0), (10_000_000, 0), (10_000_000, 10_000_000), (0, 10_000_000)]);
    let tol = ArcTol::new(1_000, Side::Outside);
    let pad = Circle::new(Point::new(3_000_000, 3_000_000), 500_000).to_ring(tol)?;
    let obstacle = offset(&pad, 200_000, Join::Round, tol)?;

    let fill = Boolean::new()
        .subject(&zone, FillRule::NonZero)
        .clip(&obstacle, FillRule::NonZero)
        .op(Op::Difference)
        .execute()?;

    let gerber_regions: Vec<Ring> = fill.iter().map(fracture).collect::<Result<_>>()?;
    assert_eq!(gerber_regions.len(), 1);
    assert!(distance_less_than(&Path::from([(0, 2_500_000), (10_000_000, 2_500_000)]), &pad, 200_000));
    Ok(())
}
```

## Performance

Single thread unless noted, Apple Silicon laptop, release build (`cargo bench`,
`examples/perf.rs`, and the Clipper2 comparison in `oracle/`):

| Workload | polyclip | with `rayon` | Clipper2 (C++) |
|---|---|---|---|
| zone 100 mm × 100 mm − 5 000 inflated obstacles | ~25 ms | ~8 ms | ~41 ms |
| union of 50 000 heavily overlapping 64-vertex circles | ~2.8 s | ~0.4 s | ~157 s |
| offset of a 10 000-vertex polygon (round joins), convex | ~1.2 ms | | ~0.9 ms |
| offset of a 10 000-vertex wavy star | ~21 ms | | |
| `distance_less_than`, 64-vertex polygons (average incl. bbox rejection) | ~80 ns | | |
| fracture of the zone result (3 672 holes) | ~10 ms | | |
| triangulation of the zone result | ~33 ms | | |
| union of 50 000 stacked slots / diagonal slots (long dense parallel edges) | ~24 / ~95 ms | | |
| cadlab GND pour (1 polygon, 2 153 holes, 242 k vertices, `testdata/`): `opening` by 100 µm | ~530 ms | ~93 ms | |
| … the opened pour ∪ 40 small rectangles (thermal spokes) | ~9 ms | ~5 ms | |
| … `offset` by −100 µm | ~76 ms | ~20 ms | |

The noder adapts its spatial index to the data (uniform grid or k-d tree, sweeps along
the direction in which the segments are thinnest, including the dominant segment
direction), and distance queries use direction-adaptive sweeps and bounding-volume
hierarchies, so long, dense or parallel edges at any angle do not degrade into quadratic
behaviour. Booleans split their input into clusters of rings that cannot interact (no edge
boxes within a unit of each other) and compute each on its own, in parallel with `rayon`;
a lone simple ring is only located and passed through, so a small operand meeting a few
rings of a large one costs little more than reading it. Within one large cluster, every
phase runs in parallel with `rayon`: noding by leaf and segment ranges, fragment merging by
`x` ranges, and the sweep in vertical bands joined at their common sides. The output is
identical to computing everything at once, on any number of threads. Offsets build their raw curves in parallel with `rayon`, and a
single raw curve that is already simple (a convex polygon grown, say) is passed through
without noding.

## Limitations

- The spec's indicative target for the union of 50 000 heavily overlapping circles
  (< 300 ms) is met only roughly, and only with `rayon` on many cores (~0.4 s on 16; ~2.8 s
  single-threaded; Clipper2 needs minutes on the same input).
  Realistic zone fills, offsets and distance queries meet their targets.
- `curved_boolean` returns real arcs (source centre and radius) with a documented
  deviation bound and side guarantee, but it is built on approximation plus
  reconstruction, not an exact line/arc arrangement: near tangencies and crowded spots
  short polyline pieces remain.

## Feature flags

- `serde`: `Serialize`/`Deserialize` for all data types.
- `rayon`: parallelize the heavy phases of boolean operations (noding, merging and the
  sweep within one cluster), and independent clusters of rings. Output is identical for any
  number of threads.

## Verification

- Property tests (`proptest`) for boolean identities, validity and canonical form,
  idempotence, offsets, distances, fracture, simplification, triangulation and
  decomposition; clustered booleans are checked bit for bit against the one-piece
  computation, and the parallel boolean engine against the original sequential one.
- Differential testing against Clipper2 (as an oracle only) in the separate `oracle/`
  crate.
- Fuzzing of every public operation with `cargo-fuzz` (`fuzz/`), run weekly in CI.
- Criterion benchmarks (`benches/`).

## MSRV

Rust 1.89.

## License

MIT
