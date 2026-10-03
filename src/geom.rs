//! Core geometric data types.

use core::ops::{Deref, DerefMut};

/// Largest supported absolute coordinate value: `2^40`.
///
/// Every coordinate given to (or produced by) an operation of this crate must satisfy
/// `|x| <= MAX_COORD && |y| <= MAX_COORD`. Operations return
/// [`Error::CoordinateOutOfRange`](crate::Error::CoordinateOutOfRange) for inputs outside
/// this range instead of risking a wrong answer. In nanometers this is about 1.1 km.
///
/// The bound guarantees that all exact predicates fit in `i128` arithmetic (or the crate's
/// internal 256-bit arithmetic for squared distances).
pub const MAX_COORD: i64 = 1 << 40;

/// A point with integer coordinates.
///
/// The derived ordering is lexicographic: by `x`, then by `y`. This is the order used for
/// canonical ring starts and for every sweep in the crate.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Point {
    /// X coordinate.
    pub x: i64,
    /// Y coordinate (the crate uses the Y-up convention: counter-clockwise means positive area).
    pub y: i64,
}

impl Point {
    /// Creates a point.
    #[inline]
    pub const fn new(x: i64, y: i64) -> Self {
        Point { x, y }
    }

    /// Returns `true` when both coordinates are within [`MAX_COORD`].
    #[inline]
    pub const fn in_range(self) -> bool {
        self.x >= -MAX_COORD && self.x <= MAX_COORD && self.y >= -MAX_COORD && self.y <= MAX_COORD
    }
}

impl From<(i64, i64)> for Point {
    #[inline]
    fn from((x, y): (i64, i64)) -> Self {
        Point { x, y }
    }
}

/// A point with floating-point coordinates, used for results that are not representable on
/// the integer grid (centroids, closest points).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PointF {
    /// X coordinate.
    pub x: f64,
    /// Y coordinate.
    pub y: f64,
}

impl PointF {
    /// Creates a point.
    #[inline]
    pub const fn new(x: f64, y: f64) -> Self {
        PointF { x, y }
    }
}

impl From<Point> for PointF {
    #[inline]
    fn from(p: Point) -> Self {
        PointF {
            x: p.x as f64,
            y: p.y as f64,
        }
    }
}

/// An axis-aligned rectangle with inclusive bounds.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Rect {
    /// Lower-left corner (inclusive).
    pub min: Point,
    /// Upper-right corner (inclusive).
    pub max: Point,
}

impl Rect {
    /// Creates a rectangle from two corners in any order.
    pub fn new(a: Point, b: Point) -> Self {
        Rect {
            min: Point::new(a.x.min(b.x), a.y.min(b.y)),
            max: Point::new(a.x.max(b.x), a.y.max(b.y)),
        }
    }

    /// Bounding box of a point set, or `None` when it is empty.
    pub fn of_points<'a, I: IntoIterator<Item = &'a Point>>(pts: I) -> Option<Rect> {
        let mut it = pts.into_iter();
        let first = *it.next()?;
        let mut r = Rect {
            min: first,
            max: first,
        };
        for p in it {
            r.add_point(*p);
        }
        Some(r)
    }

    /// Grows the rectangle to include `p`.
    #[inline]
    pub fn add_point(&mut self, p: Point) {
        self.min.x = self.min.x.min(p.x);
        self.min.y = self.min.y.min(p.y);
        self.max.x = self.max.x.max(p.x);
        self.max.y = self.max.y.max(p.y);
    }

    /// Smallest rectangle containing both.
    #[inline]
    pub fn union(&self, o: &Rect) -> Rect {
        Rect {
            min: Point::new(self.min.x.min(o.min.x), self.min.y.min(o.min.y)),
            max: Point::new(self.max.x.max(o.max.x), self.max.y.max(o.max.y)),
        }
    }

    /// `true` when the closed rectangles share at least one point.
    #[inline]
    pub fn intersects(&self, o: &Rect) -> bool {
        self.min.x <= o.max.x
            && o.min.x <= self.max.x
            && self.min.y <= o.max.y
            && o.min.y <= self.max.y
    }

    /// `true` when `p` lies in the closed rectangle.
    #[inline]
    pub fn contains_point(&self, p: Point) -> bool {
        p.x >= self.min.x && p.x <= self.max.x && p.y >= self.min.y && p.y <= self.max.y
    }

    /// `true` when `o` lies entirely inside this closed rectangle.
    #[inline]
    pub fn contains_rect(&self, o: &Rect) -> bool {
        self.contains_point(o.min) && self.contains_point(o.max)
    }

    /// The rectangle grown by `d` on every side (saturating).
    #[inline]
    pub fn expand(&self, d: i64) -> Rect {
        Rect {
            min: Point::new(self.min.x.saturating_sub(d), self.min.y.saturating_sub(d)),
            max: Point::new(self.max.x.saturating_add(d), self.max.y.saturating_add(d)),
        }
    }

    /// Width (`max.x - min.x`).
    #[inline]
    pub fn width(&self) -> i64 {
        self.max.x - self.min.x
    }

    /// Height (`max.y - min.y`).
    #[inline]
    pub fn height(&self) -> i64 {
        self.max.y - self.min.y
    }
}

macro_rules! point_vec_newtype {
    ($(#[$m:meta])* $name:ident) => {
        $(#[$m])*
        #[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
        #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
        #[cfg_attr(feature = "serde", serde(transparent))]
        pub struct $name(pub Vec<Point>);

        impl $name {
            /// Creates an empty instance.
            #[inline]
            pub const fn new() -> Self {
                $name(Vec::new())
            }

            /// Bounding box, or `None` when empty.
            #[inline]
            pub fn bbox(&self) -> Option<Rect> {
                Rect::of_points(self.0.iter())
            }

            /// Consumes `self`, returning the vertices.
            #[inline]
            pub fn into_inner(self) -> Vec<Point> {
                self.0
            }
        }

        impl Deref for $name {
            type Target = Vec<Point>;
            #[inline]
            fn deref(&self) -> &Vec<Point> {
                &self.0
            }
        }

        impl DerefMut for $name {
            #[inline]
            fn deref_mut(&mut self) -> &mut Vec<Point> {
                &mut self.0
            }
        }

        impl From<Vec<Point>> for $name {
            #[inline]
            fn from(v: Vec<Point>) -> Self {
                $name(v)
            }
        }

        impl From<&[Point]> for $name {
            #[inline]
            fn from(v: &[Point]) -> Self {
                $name(v.to_vec())
            }
        }

        impl From<&[(i64, i64)]> for $name {
            fn from(v: &[(i64, i64)]) -> Self {
                $name(v.iter().map(|&p| Point::from(p)).collect())
            }
        }

        impl<const N: usize> From<[(i64, i64); N]> for $name {
            fn from(v: [(i64, i64); N]) -> Self {
                $name(v.iter().map(|&p| Point::from(p)).collect())
            }
        }

        impl FromIterator<Point> for $name {
            fn from_iter<I: IntoIterator<Item = Point>>(it: I) -> Self {
                $name(it.into_iter().collect())
            }
        }

        impl AsRef<[Point]> for $name {
            #[inline]
            fn as_ref(&self) -> &[Point] {
                &self.0
            }
        }
    };
}

point_vec_newtype!(
    /// An open polyline: consecutive vertices are joined, the last is not joined to the first.
    Path
);

point_vec_newtype!(
    /// A closed ring: consecutive vertices are joined and an implicit edge joins the last
    /// vertex to the first. The first vertex is **not** repeated at the end.
    ///
    /// Canonical rings produced by this crate are simple (no repeated vertex, no crossing),
    /// have no zero-length or (by default) collinear-redundant edges, and start at their
    /// lexicographically smallest vertex. Outer rings are counter-clockwise, holes clockwise
    /// (Y-up convention).
    Ring
);

impl Ring {
    /// Twice the signed area, exact. Positive for counter-clockwise rings (Y-up).
    pub fn signed_area2(&self) -> i128 {
        crate::query::ring_area2(&self.0)
    }

    /// `true` when the ring is counter-clockwise (positive signed area).
    pub fn is_ccw(&self) -> bool {
        self.signed_area2() > 0
    }

    /// Reverses the orientation in place, keeping the first vertex first.
    pub fn reverse_orientation(&mut self) {
        if self.0.len() > 1 {
            self.0[1..].reverse();
        }
    }

    /// Iterates over the edges `(p[i], p[i+1])`, including the closing edge.
    pub fn edges(&self) -> impl Iterator<Item = (Point, Point)> + '_ {
        let n = self.0.len();
        (0..n).map(move |i| (self.0[i], self.0[(i + 1) % n]))
    }
}

impl From<Ring> for Path {
    /// Converts a ring into an open path that repeats the first vertex at the end (so the
    /// path traces the full closed curve).
    fn from(r: Ring) -> Path {
        let mut v = r.0;
        if let Some(&f) = v.first() {
            v.push(f);
        }
        Path(v)
    }
}

/// A polygon: one outer ring and zero or more holes.
///
/// Canonical polygons have a counter-clockwise outer ring and clockwise holes, each hole
/// strictly inside the outer ring except for isolated touching points. Holes are sorted by
/// their vertex sequence.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Polygon {
    /// Outer boundary.
    pub outer: Ring,
    /// Holes.
    pub holes: Vec<Ring>,
}

impl Polygon {
    /// Creates a polygon.
    pub fn new(outer: impl Into<Ring>, holes: Vec<Ring>) -> Self {
        Polygon {
            outer: outer.into(),
            holes,
        }
    }

    /// Iterates over the outer ring and then the holes.
    pub fn rings(&self) -> impl Iterator<Item = &Ring> {
        core::iter::once(&self.outer).chain(self.holes.iter())
    }

    /// Bounding box of the outer ring.
    pub fn bbox(&self) -> Option<Rect> {
        self.outer.bbox()
    }

    /// Twice the signed area (outer minus holes when canonical), exact.
    pub fn signed_area2(&self) -> i128 {
        self.rings().map(|r| r.signed_area2()).sum()
    }

    /// Total number of vertices.
    pub fn vertex_count(&self) -> usize {
        self.rings().map(|r| r.len()).sum()
    }
}

impl From<Ring> for Polygon {
    fn from(r: Ring) -> Self {
        Polygon {
            outer: r,
            holes: Vec::new(),
        }
    }
}

/// A set of polygons with disjoint interiors (a "multipolygon").
pub type PolygonSet = Vec<Polygon>;

/// A ring whose edges carry user tags (see [the provenance section](crate#vertex-provenance-z-tags)).
///
/// `tags[i]` is the tag of the edge from `points[i]` to `points[(i + 1) % n]`. A vertex
/// `points[i]` therefore sits between the edges tagged `tags[i - 1]` and `tags[i]`; when the
/// two differ, the vertex is where two source edges meet (often a new intersection vertex).
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TaggedRing {
    /// Vertices, as in [`Ring`].
    pub points: Vec<Point>,
    /// One tag per edge; same length as `points`.
    pub tags: Vec<u64>,
}

impl TaggedRing {
    /// Builds a tagged ring giving every edge the same tag.
    pub fn uniform(ring: impl Into<Ring>, tag: u64) -> Self {
        let r = ring.into();
        let n = r.len();
        TaggedRing {
            points: r.0,
            tags: vec![tag; n],
        }
    }

    /// Drops the tags.
    pub fn into_ring(self) -> Ring {
        Ring(self.points)
    }
}

/// An open path whose edges carry user tags. `tags[i]` belongs to the edge
/// `points[i] -> points[i + 1]`, so `tags.len() == points.len() - 1` (or 0 when empty).
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TaggedPath {
    /// Vertices.
    pub points: Vec<Point>,
    /// One tag per edge.
    pub tags: Vec<u64>,
}

impl TaggedPath {
    /// Builds a tagged path giving every edge the same tag.
    pub fn uniform(path: impl Into<Path>, tag: u64) -> Self {
        let p = path.into();
        let n = p.len().saturating_sub(1);
        TaggedPath {
            points: p.0,
            tags: vec![tag; n],
        }
    }

    /// Drops the tags.
    pub fn into_path(self) -> Path {
        Path(self.points)
    }
}

/// A polygon with tagged rings.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TaggedPolygon {
    /// Outer boundary.
    pub outer: TaggedRing,
    /// Holes.
    pub holes: Vec<TaggedRing>,
}

impl TaggedPolygon {
    /// Drops the tags.
    pub fn into_polygon(self) -> Polygon {
        Polygon {
            outer: self.outer.into_ring(),
            holes: self.holes.into_iter().map(TaggedRing::into_ring).collect(),
        }
    }
}

/// One node of a [`PolyTree`]: an outer ring or a hole.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PolyNode {
    /// The ring (counter-clockwise for outer rings, clockwise for holes).
    pub ring: Ring,
    /// Edge tags, one per edge of `ring` (all zero when the input carried no tags).
    pub tags: Vec<u64>,
    /// `true` for holes.
    pub is_hole: bool,
    /// Index of the parent node, if any. Holes always have an outer parent; outer rings
    /// either are roots or sit inside a hole (islands).
    pub parent: Option<usize>,
    /// Indices of child nodes, in canonical order.
    pub children: Vec<usize>,
}

/// Full nesting of a polygon set: outer → holes → islands inside holes → ...
///
/// Nodes live in an arena (`nodes`); `roots` lists top-level outer rings. All lists are in
/// canonical order (by ring vertex sequence), so equal regions give equal trees.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PolyTree {
    /// Node arena.
    pub nodes: Vec<PolyNode>,
    /// Indices of the top-level outer rings.
    pub roots: Vec<usize>,
}

impl PolyTree {
    /// `true` when the tree has no rings.
    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }

    /// Iterates over every polygon (each outer ring with its direct holes), in depth-first
    /// canonical order. Islands inside holes are yielded as separate polygons.
    pub fn polygons(&self) -> impl Iterator<Item = Polygon> + '_ {
        self.outer_indices().map(move |i| self.polygon_at(i))
    }

    /// Like [`polygons`](Self::polygons) but keeping edge tags.
    pub fn tagged_polygons(&self) -> impl Iterator<Item = TaggedPolygon> + '_ {
        self.outer_indices().map(move |i| {
            let n = &self.nodes[i];
            TaggedPolygon {
                outer: TaggedRing {
                    points: n.ring.0.clone(),
                    tags: n.tags.clone(),
                },
                holes: n
                    .children
                    .iter()
                    .map(|&h| TaggedRing {
                        points: self.nodes[h].ring.0.clone(),
                        tags: self.nodes[h].tags.clone(),
                    })
                    .collect(),
            }
        })
    }

    /// The polygon whose outer ring is node `i` (which must not be a hole).
    pub fn polygon_at(&self, i: usize) -> Polygon {
        let n = &self.nodes[i];
        Polygon {
            outer: n.ring.clone(),
            holes: n
                .children
                .iter()
                .map(|&h| self.nodes[h].ring.clone())
                .collect(),
        }
    }

    /// Indices of outer-ring nodes in depth-first canonical order.
    pub fn outer_indices(&self) -> impl Iterator<Item = usize> + '_ {
        let mut stack: Vec<usize> = self.roots.iter().rev().copied().collect();
        core::iter::from_fn(move || {
            while let Some(i) = stack.pop() {
                let n = &self.nodes[i];
                // Push grandchildren (islands) so they come out after this polygon.
                for &h in n.children.iter().rev() {
                    for &isl in self.nodes[h].children.iter().rev() {
                        stack.push(isl);
                    }
                }
                if !n.is_hole {
                    return Some(i);
                }
            }
            None
        })
    }

    /// Flattens the tree into a canonical [`PolygonSet`] (sorted by outer ring).
    pub fn to_polygon_set(&self) -> PolygonSet {
        let mut v: PolygonSet = self.polygons().collect();
        v.sort_by(|a, b| a.outer.0.cmp(&b.outer.0));
        v
    }

    /// Flattens the tree into canonical tagged polygons (sorted by outer ring).
    pub fn to_tagged_polygons(&self) -> Vec<TaggedPolygon> {
        let mut v: Vec<TaggedPolygon> = self.tagged_polygons().collect();
        v.sort_by(|a, b| a.outer.points.cmp(&b.outer.points));
        v
    }

    /// Twice the total signed area (exact).
    pub fn signed_area2(&self) -> i128 {
        self.nodes.iter().map(|n| n.ring.signed_area2()).sum()
    }
}
