//! Error type.

use crate::geom::Point;
use core::fmt;

/// Errors returned by fallible operations.
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum Error {
    /// An input (or computed) coordinate lies outside `±`[`MAX_COORD`](crate::MAX_COORD).
    CoordinateOutOfRange(Point),
    /// A parameter is invalid (negative tolerance, zero radius where one is required, ...).
    InvalidParameter(&'static str),
    /// The input is too large for the operation (more than `u32::MAX` edges).
    TooLarge,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::CoordinateOutOfRange(p) => {
                write!(
                    f,
                    "coordinate ({}, {}) outside the supported range ±2^40",
                    p.x, p.y
                )
            }
            Error::InvalidParameter(s) => write!(f, "invalid parameter: {s}"),
            Error::TooLarge => write!(f, "input too large"),
        }
    }
}

impl std::error::Error for Error {}

/// Result alias.
pub type Result<T> = core::result::Result<T, Error>;

#[inline]
pub(crate) fn check_point(p: Point) -> Result<()> {
    if p.in_range() {
        Ok(())
    } else {
        Err(Error::CoordinateOutOfRange(p))
    }
}
