//! Exact integer 2D polygon geometry.
//!
//! `polyclip` provides robust polygon operations on integer coordinates: boolean operations
//! (union, intersection, difference, xor) with arbitrary fill rules, N-ary unions, clipping
//! of open paths, offsetting (polygons and open paths, all join and cap types), arc and
//! circle approximation with a selectable error side, exact distance queries, point
//! location, containment and intersection predicates, validity checks, fracturing into
//! hole-free outlines, topology-preserving simplification, convex hulls and Minkowski sums.
//!
//! It was written for the CAD tool *cadlab* (copper zone fills, DRC, silkscreen clipping,
//! pad shapes, Gerber output) but has no domain-specific code.
//!
//! ```
//! use polyclip::*;
//!
//! // A 10 mm square zone (nanometer units) minus two round obstacles grown by a clearance.
//! let zone = Ring::from([(0, 0), (10_000_000, 0), (10_000_000, 10_000_000), (0, 10_000_000)]);
//! let tol = ArcTol::new(1_000, Side::Outside);
//! let pads: Vec<Ring> = [(3_000_000, 3_000_000), (7_000_000, 6_000_000)]
//!     .iter()
//!     .map(|&(x, y)| Circle::new(Point::new(x, y), 500_000).to_ring(tol).unwrap())
//!     .collect();
//! let obstacles = offset(&pads, 200_000, Join::Round, tol).unwrap();
//!
//! let fill: PolyTree = Boolean::new()
//!     .subject(&zone, FillRule::NonZero)
//!     .clip(&obstacles, FillRule::NonZero)
//!     .op(Op::Difference)
//!     .execute_tree()
//!     .unwrap();
//! assert_eq!(fill.polygons().count(), 1);
//!
//! // Minimum-width enforcement, then Gerber-ready hole-free regions.
//! let fill = opening(&fill, 100_000, ArcTol::new(1_000, Side::Inside)).unwrap();
//! let regions: Vec<Ring> = fill.iter().map(|p| fracture(p).unwrap()).collect();
//! assert_eq!(regions.len(), 1);
//!
//! // DRC: is a track closer than the clearance to a pad?
//! let track = Path::from([(0, 2_500_000), (10_000_000, 2_500_000)]);
//! assert!(distance_less_than(&track, &pads[0], 200_000));
//! ```
//!
//! # Number model
//!
//! * Coordinates are `i64` ([`Point`]); every input and output coordinate must lie within
//!   `±`[`MAX_COORD`] (`2^40`, about 1.1 km in nanometers). Out-of-range input is reported
//!   as [`Error::CoordinateOutOfRange`], never silently mishandled.
//! * All topological decisions (orientation, intersection, point location, distance
//!   comparisons) are computed exactly with `i128` arithmetic (and a 384-bit product for
//!   squared distances). Floating point is used only to construct new vertices (arcs,
//!   offsets) and as a filter in front of exact fallbacks.
//! * New vertices are integers. Boolean operations use *snap rounding*: every crossing is
//!   rounded to the nearest integer point (`floor(v + 1/2)` per coordinate) and every edge
//!   passing through the same unit pixel is routed through that point. Edges that are not
//!   moved by any rounding stay exactly where they were (so valid input passes through
//!   unchanged). No vertex moves by more than `sqrt(2)/2`, and the output is always valid:
//!   simple rings, no crossing edges, correct nesting.
//! * Arc and offset vertices are computed in `f64` (with the pure-Rust `libm`, identical on
//!   every platform) and rounded to the nearest integer point.
//!
//! # Canonical output
//!
//! Every polygon-producing operation returns canonical data, so results are deterministic
//! and comparable with `==`:
//!
//! * outer rings counter-clockwise, holes clockwise (Y-up: counter-clockwise means positive
//!   [signed area](Ring::signed_area2));
//! * no repeated consecutive vertices, no zero-length edges, collinear vertices removed
//!   (configurable on [`Boolean`]) except where needed to keep tags or shared vertices;
//! * each ring starts at its lexicographically smallest vertex (min `x`, then min `y`);
//! * holes sorted by vertex sequence, polygons sorted by their outer ring's vertex sequence;
//! * regions touching at a single point are separate rings (rings never share an edge);
//! * bit-identical across platforms and runs.
//!
//! [`check_canonical`] verifies all of the above; [`validate`] checks validity of arbitrary
//! input with a reason.
//!
//! # Fill rules and orientation
//!
//! Input rings may have any orientation, overlap and self-intersect. Winding numbers count
//! counter-clockwise rings positively; [`FillRule`] (even-odd, non-zero, positive, negative)
//! decides what is inside, separately for the subject and the clip operand.
//!
//! # Vertex provenance (Z-tags)
//!
//! Every input edge can carry a `u64` tag ([`TaggedRing`], [`TaggedPath`],
//! [`TaggedPolygon`]); untagged edges carry `0`. Output edges keep the tag of the input
//! edge they lie on, through booleans ([`Boolean::execute_tagged`],
//! [`Boolean::execute_tree`]), open path clipping ([`clip_paths`]) and offsets
//! ([`offset_tagged`]). A vertex sits between two tagged edges, so a vertex created at an
//! intersection records both source tags. Uses:
//!
//! * reconstruct arcs after booleans and offsets: tag the edges of each approximated arc
//!   with an arc id ([`Shape::to_tagged`], [`offset_shape_tagged`]); consecutive output
//!   edges with the same id form one arc (emit it as a single Gerber arc);
//! * explain results: "this fill edge comes from the clearance around U3 pad 7".
//!
//! # Robustness
//!
//! No operation panics on any input: degenerate rings, collinear or duplicate points,
//! coincident edges, zero area and huge vertex counts give a well-defined result or an
//! error. The crate has no `unsafe` code and no global state; all types are `Send + Sync`.
//! The [`Boolean`] engine can be reused to avoid reallocations.
//!
//! # Features
//!
//! * `serde`: `Serialize`/`Deserialize` for all data types.
//! * `rayon`: parallelize the heavy phases of booleans (pair search, leaf processing, large
//!   sorts). Results are identical to the sequential build, for any number of threads.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

mod arc;
mod arrangement;
mod assemble;
mod boolean;
mod decompose;
mod distance;
mod error;
mod fracture;
mod geom;
mod hull;
mod node;
mod offset;
mod par;
pub mod predicates;
mod query;
mod simplify;
mod sweep;
mod triangulate;
mod validate;
mod wide;

pub use arc::{ArcTol, Circle, Contour, Curve, MAX_ARC_VERTICES, Shape, Side};
pub use boolean::{
    Boolean, ClippedPaths, FillRule, Op, PathSource, RingSource, VertexVisitor, boolean,
    clip_paths, union_all,
};
pub use decompose::{Trapezoid, trapezoids};
pub use distance::{Closest, SqDist, distance, distance_less_than, distance_sq};
pub use error::{Error, Result};
pub use fracture::{fracture, fracture_set};
pub use geom::*;
pub use hull::{convex_hull, convex_hull_of, minkowski_sum};
pub use offset::{
    EndCap, Join, closing, offset, offset_paths, offset_paths_tagged, offset_paths_tree,
    offset_shape, offset_shape_tagged, offset_tagged, offset_tree, opening,
};
pub use query::{
    Geometry, Location, Segment, area2, centroid, contains, intersects, locate, locate_in_polygon,
    locate_in_ring, ring_area2, ring_winding,
};
pub use simplify::{simplify_path, simplify_polygon, simplify_polygons};
pub use triangulate::{Triangulation, triangulate, triangulate_delaunay, triangulate_set};
pub use validate::{RingId, ValidityError, check_canonical, validate, validate_set};
