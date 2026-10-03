# polyclip-oracle

Differential testing of `polyclip` against [Clipper2] (Angus Johnson, Boost Software
License), the reference polygon clipping library, as required by the spec
(POLYGON_LIB.md §6). Clipper2 is an **oracle only**: this crate is excluded from the
root workspace and Clipper2 is never a dependency of the polyclip library.

Clipper2 2.0.1 is used through its C API, [`clipper2c-sys`], which vendors and compiles
the original C++ sources with `cc` (a C++17 compiler is required).

[Clipper2]: https://github.com/AngusJohnson/Clipper2
[`clipper2c-sys`]: https://crates.io/crates/clipper2c-sys

## Running

```bash
cd oracle
cargo test --release                          # differential tests (< 1 s)
cargo test --release -- --nocapture           # ... with per-configuration statistics
cargo run --release --example timing          # speed comparison (union of 10 000 circles)
cargo run --release --example timing -- 50000 # full spec size (Clipper2 takes minutes)
```

CI: `.github/workflows/oracle.yml` runs `cargo test --release` here on every push to
`master` and on pull requests.

## How results are compared

Both libraries produce integer vertices, so results are never bit-identical in general:
polyclip snap-rounds (every point moves by less than one unit, edges are routed through
the unit pixels of nearby vertices), Clipper2 rounds intersection points. Results are
compared by

* the **total area**, and
* the **area of the symmetric difference** `area(xor(polyclip, clipper2))`, computed
  exactly with polyclip itself (Clipper2's output, outer rings positive / holes negative,
  is read with the non-zero rule),

both within a tolerance proportional to the perimeter:

| Comparison              | Tolerance (area units)                                     |
|-------------------------|------------------------------------------------------------|
| booleans, general input | `2 × perimeter(inputs) + 4` (rounding band on both sides)  |
| booleans, rectangles    | `0`: every intersection is integral, results must be equal |
| offsets / strokes       | `(arc tolerance + 2) × perimeter(both results) + 4`: the libraries flatten arcs with different vertex counts and phases (polyclip uses `Side::Inside`, i.e. chords, like Clipper2) |

Each test reports, per configuration, the largest symmetric difference and the largest
fraction of the tolerance used, and fails listing every case beyond tolerance with its
full input.

## What is covered (`tests/differential.rs`)

| Test                       | Inputs                                                                                     |
|----------------------------|--------------------------------------------------------------------------------------------|
| `boolean_random_rings`     | 1–3 random self-intersecting rings per operand (3–12 vertices) at coordinate ranges ±8 (heavy degeneracy), ±1000, ±10⁶ and ±2⁴⁰; every op × every fill rule |
| `boolean_rectangles_exact` | 1–8 rectangles per operand, either orientation, ranges 10 / 1000 / 2³⁰; must match exactly |
| `boolean_cad_like`         | 5–30 rectangles, circles (16/32/64 vertices) and tracks (stadiums, mostly 0/45/90°) per operand in ±10 mm (nm units); every op, non-zero and even-odd |
| `offset_random_polygons`   | normalized random polygons, both delta signs, round / square / bevel joins                |
| `offset_cad_like`          | unions of CAD-like shapes, both delta signs (up to ±0.5 mm), round / square / bevel / miter (limit 2) joins |
| `offset_paths_strokes`     | random and rectilinear polylines, every end cap × every join                              |
| `harness_conventions`      | the wrapper's orientation, fill-rule and stroke-width conventions agree with polyclip's   |

### Known, by-design differences (excluded from comparison)

* **Over-long miters.** When a miter would exceed `limit × delta`, polyclip cuts it at
  that distance (perpendicular to the bisector); Clipper2 falls back to a square join at
  distance `delta`. Miters are therefore compared only where every corner turns by at
  most 120° (miter length ≤ 2 δ), e.g. CAD shapes and rectilinear paths.
* **`EndCap::Joined` on self-crossing loops.** Clipper2 strokes a closed path by
  offsetting it outwards and inwards and filling with the positive rule, which drops
  parts of the band around a self-crossing loop (e.g. a 321.7·10⁶ nm² stroke comes out as
  230.7·10⁶); polyclip returns the full stroke, equal to the stroke of the explicitly
  closed open path. Joined strokes are compared on simple loops only.

## Timing comparison (`examples/timing.rs`)

The three spec workloads, single thread, median of several runs; Clipper2 is timed
inside its `Execute` call only (C API conversions excluded), polyclip over the whole
public call. Each run also checks that both results agree by area.

Measured on an Apple M-series laptop (other jobs running), polyclip at commit
`ef8361d`:

| Workload                                          | polyclip | Clipper2 |
|---------------------------------------------------|---------:|---------:|
| union of 10 000 circles (64 vertices, r = 1 mm)   |   221 ms |   807 ms |
| union of 50 000 circles (spec size)               |   2.35 s |    156 s |
| 100 mm zone − 5 000 obstacles (32 vertices)       |    61 ms |    41 ms |
| offset of a 10 000-vertex polygon (+0.1 mm round) |   6.7 ms |  0.89 ms |

The default union size is 10 000 circles because Clipper2 needs minutes for 50 000
heavily overlapping circles.
