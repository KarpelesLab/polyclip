//! Exact integer 2D polygon geometry.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

mod arrangement;
mod assemble;
mod boolean;
mod error;
mod geom;
mod node;
pub mod predicates;
mod query;
mod sweep;
mod validate;

pub use boolean::{
    Boolean, ClippedPaths, FillRule, Op, PathSource, RingSource, boolean, clip_paths, union_all,
};
pub use error::{Error, Result};
pub use geom::*;
pub use query::{
    Geometry, Location, Segment, area2, centroid, contains, intersects, locate, locate_in_polygon,
    locate_in_ring, ring_area2, ring_winding,
};
pub use validate::{RingId, ValidityError, check_canonical, validate, validate_set};
