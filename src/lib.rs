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
pub mod query;
mod sweep;

pub use boolean::{
    Boolean, ClippedPaths, FillRule, Op, PathSource, RingSource, boolean, clip_paths, union_all,
};
pub use error::{Error, Result};
pub use geom::*;
