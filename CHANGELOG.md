# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.0.3](https://github.com/KarpelesLab/polyclip/compare/v0.0.2...v0.0.3) - 2026-10-03

### Fixed

- simplify_polygons no longer depends on polygon order or ring start
- triangulate polygons sharing part of an edge
- exact contour orientation for curved shapes; arc direction from the ring in arcs_from_tags
- range-check every vertex (holes included) in infallible queries

### Other

- un-ignore second adversarial review regressions
- handle many edges meeting at one vertex without quadratic pair tests
- second adversarial review

## [0.0.2](https://github.com/KarpelesLab/polyclip/compare/v0.0.1...v0.0.2) - 2026-10-03

### Added

- arcs_from_tags rebuilds arcs from tagged boolean/offset output

### Fixed

- contains() on paths ignores the path's own self-crossings
- Minkowski convex shortcut no longer accepts self-intersecting stars
- keep Nearest arc approximations within tolerance; accept nearly flat arcs
- nest rings by sweep when rings touch at vertices
- never overflow on out-of-range query input; exact centroid and trapezoid areas

### Other

- operations overview and provenance helper in crate docs
- list known limitations
- refresh benchmark numbers
- remove remaining quadratic and memory cliffs
- adversarial regression tests
- bucket the arrangement sort by x before sorting buckets
- bound the candidate window in exact distance search
- update zone fill and offset timings
- move arcs_from_tags above the test module
- compute ring nesting from regions tracked during the main sweep

## [0.0.1](https://github.com/KarpelesLab/polyclip/compare/v0.0.0...v0.0.1) - 2026-10-03

### Added

- optional rayon feature for internal parallelism
- arc/circle approximation with side selection and polygon/path offsetting
- exact queries (locate, intersects, contains, area, centroid) and validity checks
- snap-rounding noder, plane sweep and boolean operations

### Fixed

- reject i64::MIN in opening/closing, range-check shapes and single-vertex paths

### Other

- refresh performance table, document rayon and verification; test Send + Sync
- trim concave offset joins at the edge intersection when safe
- satisfy clippy type_complexity in leaf precompute
- end-to-end PCB-like zone fill pipeline
- lazy, leaf-local snap rounding and assembly without re-sorting
- adaptive spatial index for noding and blocked sweep status
- Merge branch 'worktree-agent-a9b88d6dabd25a4db'
- run the Clipper2 differential tests on push and pull requests
- add differential tests against Clipper2 in a standalone oracle crate
- check that the benchmarks compile
- add criterion benchmarks for the spec workloads
- run the fuzz targets weekly and on demand
- add cargo-fuzz targets for every public operation
- hash-based ring linking, split noder phases, float pre-filters for pixel tests
- grid-based noder, single sweep-order sort, float-filtered crossing rounding
- add CI, crates.io, docs.rs and license badges
