# Fuzz tests

Coverage-guided fuzz targets for every public operation of `polyclip`, built with
[`cargo-fuzz`] + [libFuzzer]. The spec requires that the library **never panics** on any
input and that every output is valid and canonical; these targets check both.

The parent crate's stable-Rust CI does **not** run this directory. Fuzzing requires
nightly Rust; see `.github/workflows/fuzz.yml` for the weekly-cron + manual-dispatch CI
job that does run the targets.

[`cargo-fuzz`]: https://github.com/rust-fuzz/cargo-fuzz
[libFuzzer]: https://llvm.org/docs/LibFuzzer.html

## Targets

| Target         | Operations                                                         | Checked beyond "no panic"                                                                                                                                       |
|----------------|--------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `boolean`      | `Boolean` (all ops × fill rules, `keep_collinear`), `boolean`, `union_all` | out-of-range ⇔ `Err`; output passes `check_canonical`; `execute` / `execute_tree` / `execute_tagged` / `boolean()` agree; re-normalizing output stays canonical and keeps its area (± snap slack) |
| `clip_paths`   | `clip_paths`                                                       | out-of-range ⇔ `Err`; pieces have ≥ 1 edge, one tag per edge, no zero-length edge; pieces appear iff some path has an edge; no inside piece without clip edges |
| `offset`       | `offset`, `offset_tree`, `opening`, `closing` (all joins incl. NaN / < 1 / huge miter limits, all sides, invalid tolerances, huge deltas) | out-of-range input rejected; only `CoordinateOutOfRange` / `InvalidParameter` errors; output canonical; tree form agrees; growing never loses and shrinking never gains area; `delta == 0` is normalization |
| `offset_paths` | `offset_paths`, `offset_paths_tree` (all joins × all end caps)     | negative delta / out-of-range rejected; output canonical; tree form agrees; path vertices covered (round/square caps at ends, non-bevel joins inside)          |
| `arc`          | `Circle::to_ring`, `Shape::to_polygon`, `Shape::to_tagged`          | vertices in range; circles (r ≥ 16) are canonical rings on the requested `Side`; tagged/untagged agree, tags name real contour elements; shape normalizes to a canonical set |
| `query`        | `locate`, `locate_in_ring`, `locate_in_polygon`, `ring_winding`, `ring_area2`, `intersects`, `contains`, `area2`, `centroid` on points, segments, paths, rings, polygons, sets and trees | raw (invalid) input: no panic. Canonical input: `contains(a,b) ⇒ intersects(a,b)` for non-empty `b`, `intersects` symmetric, `contains(a,a)`, point `contains`/`intersects` ⇔ `locate != Outside`, vertices `OnBoundary`, `area2` = Σ ring areas, centroid inside bbox, `contains(a,b) ⇒ area(b−a) ≈ 0` |
| `validate`     | `validate`, `validate_set`, `check_canonical`                       | out-of-range ⇒ invalid; canonical ⇒ valid ⇒ every member valid; valid sets keep their area under `union_all` (± snap slack); canonical sets normalize to canonical sets; perturbed canonical sets (moved / duplicated / removed vertex, reversed / rotated ring, duplicated polygon) |

## Input model

All targets decode their bytes through `polyclip_fuzz::Gen` (`src/lib.rs`), which biases
the geometry towards the cases that break polygon code:

- a coordinate **scale** per input (re-drawn now and then per ring): `0..=3` (heavy
  degeneracy: almost everything is collinear, touching or coincident), `±16`, `±2^20`,
  the full `±MAX_COORD = ±2^40` range, and values hugging `±MAX_COORD`;
- one input in 16 may contain coordinates **just beyond** the range (`MAX_COORD + 1`,
  `2·MAX_COORD`, `i64::MIN`, `i64::MAX`, …), which fallible operations must reject with
  `Err`, never panic on;
- **duplicate** points, **collinear** continuations of the previous edge, vertices reused
  from earlier rings, axis-aligned steps, rectangles and polyclip-approximated circles;
- arc tolerances are floored at `radius / 2^14` so a single circle never needs more than
  a few hundred vertices (millions of arc points would only make the fuzzer slow), with
  occasional invalid (`< 1`) tolerances.

Queries are infallible and documented for in-range coordinates only, so the `query`
target skips inputs with out-of-range coordinates.

## Running a single target

```bash
# One-time install:
cargo +nightly install cargo-fuzz

# 60-second run (from the repository root or from fuzz/):
cargo +nightly fuzz run boolean -- -max_total_time=60
```

By default cargo-fuzz builds with debug assertions and overflow checks, which is what we
want: an arithmetic overflow is a bug. Add `-O` to fuzz the release configuration only
(useful to keep searching while a `debug_assert!` failure is being fixed).

Interesting inputs are saved in `fuzz/corpus/<target>/`; a crash input goes to
`fuzz/artifacts/<target>/crash-<hash>`.

## Running every target (smoke)

```bash
for t in $(cargo +nightly fuzz list); do
    echo "=== $t ==="
    cargo +nightly fuzz run "$t" -- -max_total_time=30 || break
done
```

## Reproducing a crash as a readable test case

Every target prints its decoded input (as Rust `Debug`) when `POLYCLIP_FUZZ_DEBUG` is set:

```bash
POLYCLIP_FUZZ_DEBUG=1 cargo +nightly fuzz run boolean fuzz/artifacts/boolean/crash-XXXX
```

Turn the printed values into a test in `tests/fuzz_regressions.rs`.

## Adding a new target

1. Add a `[[bin]]` entry to `fuzz/Cargo.toml`.
2. Create `fuzz_targets/<name>.rs`:
   ```rust
   #![no_main]
   use libfuzzer_sys::fuzz_target;
   use polyclip_fuzz::Gen;
   fuzz_target!(|data: &[u8]| {
       let mut g = Gen::new(data);
       let rings = g.rings(4, 10);
       polyclip_fuzz::dump("rings", &rings);
       let _ = polyclip::union_all(&rings, g.fill_rule());
   });
   ```
3. Add the target name to the matrix in `.github/workflows/fuzz.yml`.
4. (Optional) Drop hand-picked seed inputs in `fuzz/corpus/<name>/`.

## Triage convention

A new crash is a real bug until proven otherwise:

1. Reproduce and decode it (see above).
2. If the library is at fault, add an `#[ignore]`d test to `tests/fuzz_regressions.rs`
   describing the bug, and fix it in `src/` (separate `fix(...)` commit, un-ignoring the
   test). Add the input as a named seed `fuzz/corpus/<target>/regression-<short-desc>`.
3. If the target's check was wrong (a property that does not actually hold, e.g. exact
   area preservation under snap rounding), fix the target and document why.
4. While a bug is open, a target may carry a narrowly scoped guard that skips the
   offending input class, with a comment naming the regression test; remove it with the
   fix.

## Additional targets

| Target | Covers |
|---|---|
| `polygon_tools` | `fracture`, `triangulate`/`triangulate_delaunay`, `simplify_polygon(s)`, `trapezoids`, `convex_hull`, `minkowski_sum` (raw input: no panic; canonical input: exact area and validity checks) |
| `distance` | `distance`, `distance_less_than` on sets, paths, segments and points (symmetry, zero distance iff intersecting, threshold agrees with the full distance, closest points realize it) |
