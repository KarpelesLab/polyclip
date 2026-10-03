//! Exact integer 2D polygon geometry.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

mod arc;
mod arrangement;
mod assemble;
mod boolean;
mod distance;
mod error;
mod fracture;
mod geom;
mod hull;
mod node;
mod offset;
pub mod predicates;
mod query;
mod simplify;
mod sweep;
mod validate;
mod wide;

pub use arc::{ArcTol, Circle, Contour, Curve, MAX_ARC_VERTICES, Shape, Side};
pub use boolean::{
    Boolean, ClippedPaths, FillRule, Op, PathSource, RingSource, boolean, clip_paths, union_all,
};
pub use distance::{Closest, SqDist, distance, distance_less_than, distance_sq};
pub use error::{Error, Result};
pub use fracture::{fracture, fracture_set};
pub use geom::*;
pub use hull::{convex_hull, convex_hull_of, minkowski_sum};
pub use offset::{
    EndCap, Join, closing, offset, offset_paths, offset_paths_tree, offset_tree, opening,
};
pub use query::{
    Geometry, Location, Segment, area2, centroid, contains, intersects, locate, locate_in_polygon,
    locate_in_ring, ring_area2, ring_winding,
};
pub use simplify::{simplify_path, simplify_polygon, simplify_polygons};
pub use validate::{RingId, ValidityError, check_canonical, validate, validate_set};
