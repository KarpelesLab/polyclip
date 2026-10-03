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

Single thread, Apple Silicon laptop, release build (see `examples/perf.rs`):

| Workload | polyclip | Clipper2 (C++) |
|---|---|---|
| zone 100 mm × 100 mm − 5 000 inflated obstacles | ~50–60 ms | ~41 ms |
| union of 50 000 heavily overlapping 64-vertex circles | ~2.5 s | ~157 s |
| offset of a 10 000-vertex polygon (round joins) | ~10 ms | — |
| `distance_less_than`, two 64-vertex polygons (average incl. bbox rejection) | ~80 ns | — |
| fracture of the zone result (3 672 holes) | ~13 ms | — |

## Feature flags

- `serde`: `Serialize`/`Deserialize` for all data types.

## MSRV

Rust 1.89.

## License

MIT
