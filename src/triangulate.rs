//! Constrained triangulation of polygons with holes.
//!
//! # Algorithm
//!
//! A single plane sweep in lexicographic `(x, y)` order (the same sheared vertical sweep
//! line as every other sweep in the crate) triangulates the polygon directly, without a
//! separate monotone decomposition pass (Hertel–Mehlhorn style). The sweep status holds
//! the boundary edges crossing the sweep line; every interval between two consecutive
//! status edges that lies inside the polygon keeps the not yet triangulated part of the
//! swept region as a *reflex chain*: a polyline from the left endpoint of the interval's
//! lower edge to the left endpoint of its upper edge whose interior vertices are all reflex
//! (or flat). When a vertex arrives at the end of the interval's upper (lower) edge, it
//! cuts off triangles from the top (bottom) of the chain as long as the chain vertex is
//! strictly convex, exactly as in the classic monotone-polygon triangulation. A vertex
//! that arrives strictly inside an interval (the leftmost vertex of a hole) connects to the
//! most recently swept vertex of the interval, which is always visible, and splits the
//! interval in two; intervals separated by a vertex where two of them meet are merged.
//!
//! Vertices are unique points: rings touching each other at shared vertices (as produced
//! by the boolean operations) are handled natively, every edge around such a vertex is
//! processed in one event. The status is a list of short sorted chunks, so updates stay
//! cheap even when tens of thousands of edges cross the sweep line; the run time is
//! `O(n log n)` for `n` vertices on typical inputs (status updates cost
//! `O(CHUNK + width / CHUNK)` each).
//!
//! All decisions use the exact orientation predicate; every emitted triangle has strictly
//! positive orientation. [`triangulate_delaunay`] then applies Lawson edge flips with an
//! exact in-circle predicate (256-bit arithmetic) to make the triangulation constrained
//! Delaunay.
//!
//! Input handling is documented on [`triangulate`].

use crate::error::{Error, Result};
use crate::geom::{Point, Polygon};
use crate::predicates::{cmp_dir_halfplane, cross, orient, sub};
use crate::query::ring_area2;
use core::cmp::Ordering;
use std::collections::VecDeque;

/// A triangulation: vertices and counter-clockwise triangles indexing them.
///
/// Produced by [`triangulate`], [`triangulate_set`] and [`triangulate_delaunay`]. The
/// layout is deterministic:
///
/// * `vertices` holds every distinct input vertex used by the triangulation exactly once,
///   sorted lexicographically (by `x`, then `y`). No new vertices are created.
/// * Each triangle is counter-clockwise with strictly positive area, and starts at its
///   smallest vertex index; the triangle list is sorted.
/// * The triangles cover the polygon exactly (their doubled areas add up to
///   [`Polygon::signed_area2`] of the canonical polygon), never overlap, and every
///   boundary edge of the polygon is an edge of exactly one triangle (unless split by a
///   vertex of another ring lying on it, which valid polygons never have).
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Triangulation {
    /// Distinct vertices, sorted lexicographically.
    pub vertices: Vec<Point>,
    /// Counter-clockwise triangles as indices into `vertices`.
    pub triangles: Vec<[u32; 3]>,
}

impl Triangulation {
    /// The three corners of triangle `i` (counter-clockwise).
    ///
    /// # Panics
    ///
    /// Panics when `i` or one of the triangle's indices is out of bounds.
    pub fn triangle(&self, i: usize) -> [Point; 3] {
        let t = self.triangles[i];
        [
            self.vertices[t[0] as usize],
            self.vertices[t[1] as usize],
            self.vertices[t[2] as usize],
        ]
    }

    /// Sum of the doubled triangle areas (exact). Equals the doubled area of the
    /// triangulated region.
    pub fn area2(&self) -> i128 {
        self.triangles
            .iter()
            .map(|t| {
                match (
                    self.vertices.get(t[0] as usize),
                    self.vertices.get(t[1] as usize),
                    self.vertices.get(t[2] as usize),
                ) {
                    (Some(&a), Some(&b), Some(&c)) => orient(a, b, c),
                    _ => 0,
                }
            })
            .sum()
    }
}

/// Triangulates a polygon with holes.
///
/// Returns counter-clockwise, non-degenerate triangles covering the polygon exactly, with
/// every boundary edge as a triangle edge (see [`Triangulation`] for the output layout).
/// No vertices are added.
///
/// The triangulation is computed by a single plane sweep (Hertel–Mehlhorn style: a reflex
/// chain per interval of the sweep status, split at hole vertices and merged where
/// intervals meet), using only the exact orientation predicate. It runs in `O(n log n)`
/// for typical inputs; a zone fill with 2000 holes and 100 000 vertices takes a few tens of
/// milliseconds.
///
/// # Input handling
///
/// * Valid polygons (as produced by the boolean operations) are fully supported,
///   including rings touching each other at shared vertices (a hole touching the outer
///   ring, holes touching each other) and collinear consecutive vertices, which are kept
///   as vertices of the triangulation.
/// * Rings may have any orientation: the outer ring is treated as counter-clockwise and
///   the holes as clockwise (each ring is reversed when needed, based on its signed area).
/// * Consecutive duplicate vertices are ignored; rings with fewer than three distinct
///   vertices or with zero signed area are skipped (they enclose no area).
/// * Identical points (in one ring or in different rings) become a single vertex.
/// * A vertex lying in the interior of another ring's edge splits that edge (the edge then
///   appears as several triangle edges).
/// * Two edges along the same segment in opposite directions cancel each other (zero-width
///   slits, a hole equal to the outer ring, polygons of a set sharing an edge).
/// * Other invalid input (crossing or overlapping edges, holes outside their outer ring,
///   overlapping holes or polygons) is detected by consistency checks of the sweep and of
///   the result (the triangle areas must add up to the enclosed area exactly) and gives
///   an error. Repair such input with [`union_all`](crate::union_all) first.
/// * No input causes a panic.
///
/// ```
/// use polyclip::{triangulate, Polygon, Ring};
/// let outer = Ring::from([(0, 0), (10, 0), (10, 10), (0, 10)]);
/// let hole = Ring::from([(3, 3), (3, 7), (7, 7), (7, 3)]);
/// let t = triangulate(&Polygon::new(outer, vec![hole])).unwrap();
/// assert_eq!(t.vertices.len(), 8);
/// assert_eq!(t.triangles.len(), 8);
/// assert_eq!(t.area2(), 2 * (100 - 16));
/// ```
///
/// # Errors
///
/// [`Error::CoordinateOutOfRange`] for coordinates beyond `±`[`MAX_COORD`](crate::MAX_COORD),
/// [`Error::TooLarge`] for more than `u32::MAX` vertices and [`Error::InvalidParameter`]
/// for invalid polygons (crossing or overlapping edges, bad nesting).
pub fn triangulate(poly: &Polygon) -> Result<Triangulation> {
    triangulate_set(core::slice::from_ref(poly))
}

/// Triangulates a set of polygons with pairwise disjoint interiors (they may touch at
/// vertices) into one triangulation sharing vertices between polygons.
///
/// Same input handling, guarantees and errors as [`triangulate`]; each polygon's first
/// ring is its outer ring. Islands inside holes of another polygon are fine.
pub fn triangulate_set(polys: &[Polygon]) -> Result<Triangulation> {
    let mesh = Mesh::build(polys)?;
    let tris = mesh.sweep(CHUNK)?;
    Ok(finish(mesh.verts, tris))
}

/// Constrained Delaunay triangulation of a polygon with holes.
///
/// Starts from [`triangulate`] and flips non-boundary edges (Lawson's algorithm) until
/// every edge is locally Delaunay: no vertex of a triangle lies strictly inside the
/// circumcircle of a triangle sharing an edge with it, unless that edge is a boundary
/// edge. This maximizes the minimum angle among triangulations with the same vertices
/// and boundary edges, so triangles are better shaped for rendering and meshing. The
/// in-circle test is exact. For co-circular vertices, the choice of diagonal is the one
/// of the initial triangulation (deterministic).
///
/// Same input handling, guarantees and errors as [`triangulate`]. Lawson's algorithm is
/// quadratic in the worst case, but typical inputs need a small multiple of `n` flips
/// (about twice the time of [`triangulate`]).
pub fn triangulate_delaunay(poly: &Polygon) -> Result<Triangulation> {
    let mesh = Mesh::build(core::slice::from_ref(poly))?;
    let mut tris = mesh.sweep(CHUNK)?;
    delaunay_flip(&mesh.verts, &mut tris);
    Ok(finish(mesh.verts, tris))
}

const INVALID: Error = Error::InvalidParameter(
    "triangulate: invalid polygon (crossing or overlapping edges, or bad ring nesting)",
);

/// A boundary edge between vertices `lo < hi`; `above` is `true` when the polygon interior
/// lies above it in the sweep status (the edge is directed `lo -> hi`).
#[derive(Clone, Copy, Debug, Default)]
struct Edge {
    lo: u32,
    hi: u32,
    above: bool,
}

/// Unique vertices and boundary edges of the input.
struct Mesh {
    verts: Vec<Point>,
    /// Edges sorted by `lo`, then bottom to top around `lo`.
    edges: Vec<Edge>,
    /// Exact doubled area enclosed by the edges.
    area2: i128,
}

impl Mesh {
    fn build(polys: &[Polygon]) -> Result<Mesh> {
        let mut total = 0usize;
        for poly in polys {
            for ring in poly.rings() {
                for &p in ring.iter() {
                    if !p.in_range() {
                        return Err(Error::CoordinateOutOfRange(p));
                    }
                }
                total += ring.len();
            }
        }
        if total > u32::MAX as usize {
            return Err(Error::TooLarge);
        }
        // Directed edges, interior on the left.
        let mut raw: Vec<(Point, Point)> = Vec::with_capacity(total);
        let mut ring_buf: Vec<Point> = Vec::new();
        for poly in polys {
            for (ri, ring) in poly.rings().enumerate() {
                ring_buf.clear();
                for &p in ring.iter() {
                    if ring_buf.last() != Some(&p) {
                        ring_buf.push(p);
                    }
                }
                while ring_buf.len() > 1 && ring_buf.first() == ring_buf.last() {
                    ring_buf.pop();
                }
                if ring_buf.len() < 3 {
                    continue;
                }
                let a = ring_area2(&ring_buf);
                if a == 0 {
                    continue;
                }
                let flip = (ri == 0) != (a > 0);
                let n = ring_buf.len();
                for i in 0..n {
                    let (p, q) = (ring_buf[i], ring_buf[(i + 1) % n]);
                    raw.push(if flip { (q, p) } else { (p, q) });
                }
            }
        }
        // Split edges at vertices lying on their interiors (exactly), so that partially
        // shared edges become identical pieces that cancel below. Proper crossings mean
        // invalid input.
        let raw: Vec<(Point, Point)> = match crate::node::node_exact(&raw) {
            Ok(frags) => frags.into_iter().map(|f| (f.a, f.b)).collect(),
            Err(_) => return Err(INVALID),
        };
        // Vertex ids (sorted unique points) of every edge start in one sort; an edge's end
        // is usually the next edge's start (rings are listed in order).
        let mut keys: Vec<(Point, u32)> = raw
            .iter()
            .enumerate()
            .map(|(k, e)| (e.0, k as u32))
            .collect();
        keys.sort_unstable();
        let mut verts: Vec<Point> = Vec::with_capacity(raw.len());
        let mut ida = vec![0u32; raw.len()];
        for &(p, k) in &keys {
            if verts.last() != Some(&p) {
                verts.push(p);
            }
            ida[k as usize] = verts.len() as u32 - 1;
        }
        drop(keys);
        let id = |p: Point| verts.binary_search(&p).unwrap_or(0) as u32;
        let und: Vec<Edge> = raw
            .iter()
            .enumerate()
            .map(|(k, &(_, b))| {
                let ib = match raw.get(k + 1) {
                    Some(n) if n.0 == b => ida[k + 1],
                    _ => id(b),
                };
                let (ia, ib) = (ida[k], ib);
                if ia < ib {
                    Edge {
                        lo: ia,
                        hi: ib,
                        above: true,
                    }
                } else {
                    Edge {
                        lo: ib,
                        hi: ia,
                        above: false,
                    }
                }
            })
            .collect();
        drop(raw);
        // Sorted by (lo, hi): counting sort by `lo`, then each bucket by `hi`.
        let und = bucket_sort(
            und,
            verts.len(),
            |e| e.lo as usize,
            |b| b.sort_unstable_by_key(|e| e.hi),
        );
        // Cancel opposite copies of the same segment; same-direction copies overlap.
        let mut edges: Vec<Edge> = Vec::with_capacity(und.len());
        let mut i = 0;
        while i < und.len() {
            let mut j = i;
            let mut net = 0i64;
            while j < und.len() && und[j].lo == und[i].lo && und[j].hi == und[i].hi {
                net += if und[j].above { 1 } else { -1 };
                j += 1;
            }
            match net {
                0 => {}
                1 | -1 => edges.push(Edge {
                    above: net > 0,
                    ..und[i]
                }),
                _ => return Err(INVALID),
            }
            i = j;
        }
        drop(und);
        let dir = |e: &Edge| sub(verts[e.hi as usize], verts[e.lo as usize]);
        // Already sorted by `lo`: order every run by direction (stable).
        for run in edges.chunk_by_mut(|a, b| a.lo == b.lo) {
            run.sort_by(|a, b| cmp_dir_halfplane(dir(a), dir(b)));
        }
        let mut area2: i128 = 0;
        for e in &edges {
            let c = cross(verts[e.lo as usize], verts[e.hi as usize]);
            area2 += if e.above { c } else { -c };
        }
        Ok(Mesh {
            verts,
            edges,
            area2,
        })
    }

    /// Runs the triangulating sweep (`chunk`: status chunk length, see [`Status`]).
    fn sweep(&self, chunk: usize) -> Result<Vec<[u32; 3]>> {
        let pts = &self.verts[..];
        let n = pts.len();
        let mut in_count = vec![0u32; n];
        for e in &self.edges {
            in_count[e.hi as usize] += 1;
        }
        let mut tris: Vec<[u32; 3]> = Vec::with_capacity(self.edges.len());
        let mut status = Status {
            chunks: Vec::new(),
            chunk: chunk.max(1),
        };
        let mut regions = Regions::default();
        let mut outs: Vec<Edge> = Vec::new();
        let mut run: Vec<Active> = Vec::new();
        let mut new_entries: Vec<Active> = Vec::new();
        let mut ei = 0usize;
        for v in 0..n as u32 {
            let vp = pts[v as usize];
            outs.clear();
            while ei < self.edges.len() && self.edges[ei].lo == v {
                outs.push(self.edges[ei]);
                ei += 1;
            }
            let k_in = in_count[v as usize] as usize;
            if outs.is_empty() && k_in == 0 {
                continue;
            }
            let seg = |e: &Active| (pts[e.lo as usize], pts[e.hi as usize]);
            let pos = status.locate(|e| {
                let (a, b) = seg(e);
                e.hi != v && orient(a, b, vp) > 0
            });
            // Edges ending at `v`, and edges passing through `v` (split there).
            run.clear();
            let mut end = pos;
            let mut ending = 0usize;
            let mut split = false;
            while let Some(&e) = status.get(end) {
                if e.hi == v {
                    ending += 1;
                } else {
                    let (a, b) = seg(&e);
                    if orient(a, b, vp) != 0 {
                        break;
                    }
                    outs.push(Edge {
                        lo: v,
                        hi: e.hi,
                        above: e.above,
                    });
                    split = true;
                }
                run.push(e);
                end = status.next(end);
            }
            if ending != k_in {
                return Err(INVALID);
            }
            let dir = |e: &Edge| sub(pts[e.hi as usize], vp);
            if split {
                outs.sort_by(|a, b| cmp_dir_halfplane(dir(a), dir(b)));
            }
            for w in outs.windows(2) {
                if cmp_dir_halfplane(dir(&w[0]), dir(&w[1])) != Ordering::Less {
                    return Err(INVALID);
                }
            }
            // Regions alternate inside / outside: check both the old and new edges.
            let below_cur = status.prev(pos);
            let below_in = below_cur
                .and_then(|c| status.get(c))
                .is_some_and(|e| e.above);
            let above_flag = status.get(end).is_none_or(|e| e.above);
            if !alternates(below_in, run.iter().map(|e| e.above), above_flag)
                || !alternates(below_in, outs.iter().map(|e| e.above), above_flag)
            {
                return Err(INVALID);
            }

            let below = if below_in {
                below_cur.and_then(|c| status.get(c)).map(|e| e.region)
            } else {
                None
            };
            // The inside region continuing above `v` after the event (regions are owned by
            // their lower edge).
            let mut upper: Option<u32> = None;
            if !run.is_empty() {
                if let Some(b) = below {
                    add_upper(regions.chain(b, v), v, pts, &mut tris);
                }
                for (i, e) in run.iter().enumerate() {
                    if !e.above {
                        continue;
                    }
                    if i + 1 < run.len() {
                        // Closed between two edges ending at `v`.
                        let ch = regions.chain(e.region, v);
                        add_upper(ch, v, pts, &mut tris);
                        if ch.len() != 2 {
                            return Err(INVALID);
                        }
                        regions.free(e.region);
                    } else {
                        add_lower(regions.chain(e.region, v), v, pts, &mut tris);
                        upper = Some(e.region);
                    }
                }
                // No edge leaves `v`: the regions below and above `v` become one.
                if outs.is_empty()
                    && let (Some(b), Some(a)) = (below, upper)
                {
                    let m = regions.merge(b, a);
                    if let Some(e) = below_cur.and_then(|c| status.get_mut(c)) {
                        e.region = m;
                    }
                }
            } else if let Some(r) = below {
                // `v` lies inside the region: connect it to the newest chain vertex.
                let u = regions.split(r);
                add_upper(regions.chain(r, v), v, pts, &mut tris);
                add_lower(regions.chain(u, v), v, pts, &mut tris);
                upper = Some(u);
            }
            new_entries.clear();
            let last = outs.len().wrapping_sub(1);
            for (i, e) in outs.iter().enumerate() {
                let region = match (e.above, upper) {
                    (false, _) => u32::MAX,
                    (true, Some(u)) if i == last => u,
                    (true, _) => regions.alloc(v),
                };
                new_entries.push(Active {
                    lo: v,
                    hi: e.hi,
                    above: e.above,
                    region,
                });
            }
            status.replace(pos, run.len(), &new_entries);
        }
        if !status.is_empty() {
            return Err(INVALID);
        }
        let mut sum: i128 = 0;
        for t in &tris {
            sum += orient(pts[t[0] as usize], pts[t[1] as usize], pts[t[2] as usize]);
        }
        if sum != self.area2 {
            return Err(INVALID);
        }
        Ok(tris)
    }
}

/// `true` when `first, mid..., last` strictly alternate.
fn alternates(first: bool, mid: impl Iterator<Item = bool>, last: bool) -> bool {
    let mut prev = first;
    for f in mid {
        if f == prev {
            return false;
        }
        prev = f;
    }
    prev != last
}

/// A status entry: an edge crossing the sweep line.
#[derive(Clone, Copy, Debug)]
struct Active {
    lo: u32,
    hi: u32,
    /// Interior above (the region above this edge is inside and owned by it).
    above: bool,
    /// Region index when `above`.
    region: u32,
}

/// Default chunk length of the [`Status`] list.
const CHUNK: usize = 128;

/// The sweep status: active edges bottom to top, stored as a list of short sorted chunks
/// so that insertions and removals cost `O(chunk + width / chunk)` instead of `O(width)`.
/// No chunk is empty.
struct Status {
    chunks: Vec<Vec<Active>>,
    /// Chunks are split when longer than `2 * chunk`.
    chunk: usize,
}

/// A position in the [`Status`]: `o < chunks[c].len()`, or `c == chunks.len()` and
/// `o == 0` for the end.
#[derive(Clone, Copy, Debug)]
struct Cur {
    c: usize,
    o: usize,
}

impl Status {
    fn is_empty(&self) -> bool {
        self.chunks.is_empty()
    }

    /// First position whose entry does not satisfy `below` (which must hold for a prefix).
    fn locate(&self, mut below: impl FnMut(&Active) -> bool) -> Cur {
        let c = self
            .chunks
            .partition_point(|ch| ch.last().is_some_and(&mut below));
        match self.chunks.get(c) {
            Some(ch) => Cur {
                c,
                o: ch.partition_point(below),
            },
            None => Cur { c, o: 0 },
        }
    }

    fn get(&self, at: Cur) -> Option<&Active> {
        self.chunks.get(at.c)?.get(at.o)
    }

    fn get_mut(&mut self, at: Cur) -> Option<&mut Active> {
        self.chunks.get_mut(at.c)?.get_mut(at.o)
    }

    fn next(&self, at: Cur) -> Cur {
        match self.chunks.get(at.c) {
            Some(ch) if at.o + 1 < ch.len() => Cur {
                c: at.c,
                o: at.o + 1,
            },
            Some(_) => Cur { c: at.c + 1, o: 0 },
            None => at,
        }
    }

    fn prev(&self, at: Cur) -> Option<Cur> {
        if at.o > 0 {
            Some(Cur {
                c: at.c,
                o: at.o - 1,
            })
        } else if at.c > 0 {
            let c = at.c - 1;
            let len = self.chunks.get(c).map_or(0, Vec::len);
            Some(Cur {
                c,
                o: len.checked_sub(1)?,
            })
        } else {
            None
        }
    }

    /// Replaces the `remove` entries starting at `at` with `items`.
    fn replace(&mut self, at: Cur, remove: usize, items: &[Active]) {
        let mut left = remove;
        let mut d = at.c;
        let mut o = at.o;
        while left > 0 {
            let Some(ch) = self.chunks.get_mut(d) else {
                break;
            };
            let take = left.min(ch.len().saturating_sub(o));
            ch.drain(o..o + take);
            left -= take;
            d += 1;
            o = 0;
        }
        // Drop chunks emptied by the removal. If the chunk at `at.c` went away, `at` (with
        // `at.o == 0`) now designates the front of the following chunk, which is the same
        // position in the sequence.
        let hi = d.min(self.chunks.len());
        for i in (at.c.min(hi)..hi).rev() {
            if self.chunks[i].is_empty() {
                self.chunks.remove(i);
            }
        }
        if items.is_empty() {
            return;
        }
        let (c, o) = if let Some(ch) = self.chunks.get(at.c) {
            (at.c, at.o.min(ch.len()))
        } else if let Some(last) = self.chunks.last() {
            (self.chunks.len() - 1, last.len())
        } else {
            self.chunks.push(Vec::new());
            (0, 0)
        };
        let ch = &mut self.chunks[c];
        ch.splice(o..o, items.iter().copied());
        if ch.len() > 2 * self.chunk {
            let tail = ch.split_off(ch.len() / 2);
            self.chunks.insert(c + 1, tail);
        }
    }
}

/// A pending region: the reflex chain from the lower edge's left endpoint to the upper
/// edge's left endpoint, and its most recently swept vertex.
#[derive(Default)]
struct Region {
    chain: VecDeque<u32>,
    newest: u32,
}

/// Slab of regions with reuse of freed slots (and their allocations).
#[derive(Default)]
struct Regions {
    slots: Vec<Region>,
    free: Vec<u32>,
}

impl Regions {
    fn alloc(&mut self, v: u32) -> u32 {
        let i = match self.free.pop() {
            Some(i) => i,
            None => {
                self.slots.push(Region::default());
                (self.slots.len() - 1) as u32
            }
        };
        let r = &mut self.slots[i as usize];
        r.chain.clear();
        r.chain.push_back(v);
        r.newest = v;
        i
    }

    fn free(&mut self, i: u32) {
        self.free.push(i);
    }

    /// The chain of region `i`, to which the caller adds vertex `v` (now the newest).
    fn chain(&mut self, i: u32, v: u32) -> &mut VecDeque<u32> {
        let r = &mut self.slots[i as usize];
        r.newest = v;
        &mut r.chain
    }

    /// Splits region `r` at its newest vertex `c`: `r` keeps the chain up to `c`, the
    /// returned new region gets the chain from `c` on.
    fn split(&mut self, r: u32) -> u32 {
        let u = self.alloc(0);
        let Some((rs, us)) = two_mut(&mut self.slots, r as usize, u as usize) else {
            return u;
        };
        let newest = rs.newest;
        // Find the newest vertex scanning from both ends (cost of the smaller part).
        let len = rs.chain.len();
        let mut j = 0;
        for k in 0..len {
            let (a, b) = (k, len - 1 - k);
            if rs.chain[a] == newest {
                j = a;
                break;
            }
            if rs.chain[b] == newest {
                j = b;
                break;
            }
        }
        us.chain.clear();
        if j + 1 >= len - j {
            // Upper part is the smaller one.
            us.chain.extend(rs.chain.drain(j..));
            rs.chain.push_back(newest);
        } else {
            us.chain.extend(rs.chain.drain(..=j));
            core::mem::swap(&mut rs.chain, &mut us.chain);
            us.chain.push_front(newest);
        }
        u
    }

    /// Merges region `a` (above) into region `b` (below); both chains end/start at the
    /// shared current vertex. Returns the index of the merged region.
    fn merge(&mut self, b: u32, a: u32) -> u32 {
        let Some((bs, as_)) = two_mut(&mut self.slots, b as usize, a as usize) else {
            return b;
        };
        let (keep, gone) = if bs.chain.len() >= as_.chain.len() {
            bs.chain.extend(as_.chain.iter().skip(1));
            (b, a)
        } else {
            for &x in bs.chain.iter().rev().skip(1) {
                as_.chain.push_front(x);
            }
            (a, b)
        };
        self.free(gone);
        keep
    }
}

/// Two distinct mutable elements of a slice (`None` when `i == j` or out of bounds).
fn two_mut<T>(s: &mut [T], i: usize, j: usize) -> Option<(&mut T, &mut T)> {
    if i.max(j) >= s.len() || i == j {
        return None;
    }
    Some(if i < j {
        let (l, r) = s.split_at_mut(j);
        (&mut l[i], &mut r[0])
    } else {
        let (l, r) = s.split_at_mut(i);
        (&mut r[0], &mut l[j])
    })
}

/// `v` arrives at the right end of the region's upper boundary (joined to the chain's last
/// vertex): cuts off convex chain vertices from the top.
fn add_upper(chain: &mut VecDeque<u32>, v: u32, pts: &[Point], tris: &mut Vec<[u32; 3]>) {
    let pv = pts[v as usize];
    while chain.len() >= 2 {
        let m = chain.len() - 1;
        let (a, b) = (chain[m - 1], chain[m]);
        if orient(pts[a as usize], pts[b as usize], pv) < 0 {
            tris.push([a, v, b]);
            chain.pop_back();
        } else {
            break;
        }
    }
    chain.push_back(v);
}

/// `v` arrives at the right end of the region's lower boundary (joined to the chain's first
/// vertex): cuts off convex chain vertices from the bottom.
fn add_lower(chain: &mut VecDeque<u32>, v: u32, pts: &[Point], tris: &mut Vec<[u32; 3]>) {
    let pv = pts[v as usize];
    while chain.len() >= 2 {
        let (a, b) = (chain[0], chain[1]);
        if orient(pv, pts[a as usize], pts[b as usize]) < 0 {
            tris.push([v, b, a]);
            chain.pop_front();
        } else {
            break;
        }
    }
    chain.push_front(v);
}

/// Drops unused vertices and puts the triangles in canonical order.
fn finish(verts: Vec<Point>, mut tris: Vec<[u32; 3]>) -> Triangulation {
    let mut used = vec![false; verts.len()];
    for t in &tris {
        for &i in t {
            used[i as usize] = true;
        }
    }
    let (vertices, remap) = if used.iter().all(|&u| u) {
        (verts, None)
    } else {
        let mut remap = vec![0u32; verts.len()];
        let mut out = Vec::with_capacity(verts.len());
        for (i, p) in verts.into_iter().enumerate() {
            if used[i] {
                remap[i] = out.len() as u32;
                out.push(p);
            }
        }
        (out, Some(remap))
    };
    for t in &mut tris {
        if let Some(r) = &remap {
            for i in t.iter_mut() {
                *i = r[*i as usize];
            }
        }
        // Rotate to the smallest index first.
        let k = (0..3).min_by_key(|&k| t[k]).unwrap_or(0);
        *t = [t[k], t[(k + 1) % 3], t[(k + 2) % 3]];
    }
    let n = vertices.len();
    let tris = bucket_sort(tris, n, |t| t[0] as usize, |b| b.sort_unstable());
    Triangulation {
        vertices,
        triangles: tris,
    }
}

/// Sorts `v` whose order starts with `key` (in `0..n`): a counting sort by `key`, then
/// `sort` on every bucket. Gives what a full sort by that order gives.
fn bucket_sort<T: Copy + Default>(
    v: Vec<T>,
    n: usize,
    key: impl Fn(&T) -> usize,
    sort: impl Fn(&mut [T]),
) -> Vec<T> {
    let mut start = vec![0u32; n + 1];
    for x in &v {
        start[key(x).min(n.saturating_sub(1)) + 1] += 1;
    }
    for i in 0..n {
        start[i + 1] += start[i];
    }
    let mut out = vec![T::default(); v.len()];
    let mut pos = start.clone();
    for x in v {
        let k = key(&x).min(n.saturating_sub(1));
        out[pos[k] as usize] = x;
        pos[k] += 1;
    }
    for i in 0..n {
        sort(&mut out[start[i] as usize..start[i + 1] as usize]);
    }
    out
}

const NONE: u32 = u32::MAX;

/// Lawson flips until every non-boundary edge is locally Delaunay.
fn delaunay_flip(pts: &[Point], tris: &mut [[u32; 3]]) {
    let nt = tris.len();
    // nbr[t][i]: triangle across the edge opposite corner i, NONE on the boundary.
    let mut nbr = vec![[NONE; 3]; nt];
    // Half-edges `(lo, hi, triangle, corner)` in sorted order: bucketed by `lo` (counting
    // sort), then each bucket sorted.
    let nv = tris
        .iter()
        .flatten()
        .map(|&v| v as usize + 1)
        .max()
        .unwrap_or(0);
    let mut start = vec![0u32; nv + 1];
    for tri in tris.iter() {
        for i in 0..3 {
            let (a, b) = (tri[(i + 1) % 3], tri[(i + 2) % 3]);
            start[a.min(b) as usize + 1] += 1;
        }
    }
    for v in 0..nv {
        start[v + 1] += start[v];
    }
    let mut half: Vec<(u32, u32, u32, u8)> = vec![(0, 0, 0, 0); nt * 3];
    {
        let mut pos = start.clone();
        for (t, tri) in tris.iter().enumerate() {
            for i in 0..3 {
                let (a, b) = (tri[(i + 1) % 3], tri[(i + 2) % 3]);
                let lo = a.min(b) as usize;
                half[pos[lo] as usize] = (a.min(b), a.max(b), t as u32, i as u8);
                pos[lo] += 1;
            }
        }
    }
    for v in 0..nv {
        half[start[v] as usize..start[v + 1] as usize].sort_unstable();
    }
    drop(start);
    let mut stack: Vec<(u32, u8)> = Vec::new();
    let mut i = 0;
    while i < half.len() {
        let mut j = i;
        while j < half.len() && half[j].0 == half[i].0 && half[j].1 == half[i].1 {
            j += 1;
        }
        if j - i == 2 {
            let (_, _, t, a) = half[i];
            let (_, _, u, b) = half[i + 1];
            nbr[t as usize][a as usize] = u;
            nbr[u as usize][b as usize] = t;
            stack.push((t, a));
        }
        i = j;
    }
    drop(half);
    stack.reverse();
    while let Some((t, i)) = stack.pop() {
        let (t, i) = (t as usize, i as usize);
        let u = nbr[t][i];
        if u == NONE {
            continue;
        }
        let u = u as usize;
        let Some(j) = (0..3).find(|&j| nbr[u][j] == t as u32) else {
            continue;
        };
        let (a, b, c) = (tris[t][i], tris[t][(i + 1) % 3], tris[t][(i + 2) % 3]);
        let d = tris[u][j];
        if tris[u][(j + 1) % 3] != c || tris[u][(j + 2) % 3] != b {
            continue;
        }
        let (pa, pb, pc, pd) = (
            pts[a as usize],
            pts[b as usize],
            pts[c as usize],
            pts[d as usize],
        );
        if incircle(pa, pb, pc, pd) <= 0 || orient(pa, pb, pd) <= 0 || orient(pa, pd, pc) <= 0 {
            continue;
        }
        let n_ab = nbr[t][(i + 2) % 3];
        let n_ca = nbr[t][(i + 1) % 3];
        let n_bd = nbr[u][(j + 1) % 3];
        let n_dc = nbr[u][(j + 2) % 3];
        tris[t] = [a, b, d];
        nbr[t] = [n_bd, u as u32, n_ab];
        tris[u] = [a, d, c];
        nbr[u] = [n_dc, n_ca, t as u32];
        if n_bd != NONE {
            for x in nbr[n_bd as usize].iter_mut() {
                if *x == u as u32 {
                    *x = t as u32;
                }
            }
        }
        if n_ca != NONE {
            for x in nbr[n_ca as usize].iter_mut() {
                if *x == t as u32 {
                    *x = u as u32;
                }
            }
        }
        stack.push((t as u32, 0));
        stack.push((t as u32, 2));
        stack.push((u as u32, 0));
        stack.push((u as u32, 1));
    }
}

/// Exact in-circle test: positive when `d` lies strictly inside the circumcircle of the
/// counter-clockwise triangle `abc`, zero when co-circular, negative outside.
fn incircle(a: Point, b: Point, c: Point, d: Point) -> i32 {
    let adx = (a.x - d.x) as i128;
    let ady = (a.y - d.y) as i128;
    let bdx = (b.x - d.x) as i128;
    let bdy = (b.y - d.y) as i128;
    let cdx = (c.x - d.x) as i128;
    let cdy = (c.y - d.y) as i128;
    // Floating-point filter (Shewchuk's bound for the translated determinant). The
    // differences (|v| <= 2^41) are exact in f64.
    {
        let (fadx, fady, fbdx, fbdy, fcdx, fcdy) = (
            adx as f64, ady as f64, bdx as f64, bdy as f64, cdx as f64, cdy as f64,
        );
        let alift = fadx * fadx + fady * fady;
        let blift = fbdx * fbdx + fbdy * fbdy;
        let clift = fcdx * fcdx + fcdy * fcdy;
        let bc = fbdx * fcdy - fbdy * fcdx;
        let ca = fcdx * fady - fcdy * fadx;
        let ab = fadx * fbdy - fady * fbdx;
        let det = alift * bc + blift * ca + clift * ab;
        let perm = alift * ((fbdx * fcdy).abs() + (fbdy * fcdx).abs())
            + blift * ((fcdx * fady).abs() + (fcdy * fadx).abs())
            + clift * ((fadx * fbdy).abs() + (fady * fbdx).abs());
        let bound = 4e-15 * perm;
        if det > bound {
            return 1;
        }
        if det < -bound {
            return -1;
        }
    }
    let alift = adx * adx + ady * ady;
    let blift = bdx * bdx + bdy * bdy;
    let clift = cdx * cdx + cdy * cdy;
    let bc = bdx * cdy - bdy * cdx;
    let ca = cdx * ady - cdy * adx;
    let ab = adx * bdy - ady * bdx;
    // Differences below 2^30: lifts and cross products below 2^61, so the determinant
    // (below 2^124) is exact in `i128`.
    const SMALL: i128 = 1 << 30;
    if [adx, ady, bdx, bdy, cdx, cdy]
        .iter()
        .all(|v| v.abs() < SMALL)
    {
        return (alift * bc + blift * ca + clift * ab).signum() as i32;
    }
    incircle_wide(alift, blift, clift, bc, ca, ab)
}

/// The in-circle determinant `alift * bc + blift * ca + clift * ab` in 256 bits.
fn incircle_wide(alift: i128, blift: i128, clift: i128, bc: i128, ca: i128, ab: i128) -> i32 {
    let s = I256::mul(alift, bc)
        .add(I256::mul(blift, ca))
        .add(I256::mul(clift, ab));
    s.signum()
}

/// Minimal 256-bit two's complement integer for the exact in-circle determinant.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct I256 {
    hi: u128,
    lo: u128,
}

impl I256 {
    /// Full product of two `i128` values whose magnitudes are below `2^127`.
    fn mul(a: i128, b: i128) -> I256 {
        let neg = (a < 0) != (b < 0);
        let (x, y) = (a.unsigned_abs(), b.unsigned_abs());
        const M: u128 = u64::MAX as u128;
        let (x0, x1) = (x & M, x >> 64);
        let (y0, y1) = (y & M, y >> 64);
        let p00 = x0 * y0;
        let p01 = x0 * y1;
        let p10 = x1 * y0;
        let p11 = x1 * y1;
        let mid = (p00 >> 64) + (p01 & M) + (p10 & M);
        let lo = (p00 & M) | (mid << 64);
        let hi = p11
            .wrapping_add(p01 >> 64)
            .wrapping_add(p10 >> 64)
            .wrapping_add(mid >> 64);
        let r = I256 { hi, lo };
        if neg { r.neg() } else { r }
    }

    fn neg(self) -> I256 {
        let lo = (!self.lo).wrapping_add(1);
        let hi = (!self.hi).wrapping_add(u128::from(lo == 0));
        I256 { hi, lo }
    }

    fn add(self, o: I256) -> I256 {
        let (lo, carry) = self.lo.overflowing_add(o.lo);
        let hi = self.hi.wrapping_add(o.hi).wrapping_add(u128::from(carry));
        I256 { hi, lo }
    }

    fn signum(self) -> i32 {
        if (self.hi as i128) < 0 {
            -1
        } else if self.hi == 0 && self.lo == 0 {
            0
        } else {
            1
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::Ring;
    use std::collections::HashMap;

    fn p(x: i64, y: i64) -> Point {
        Point::new(x, y)
    }

    fn ring(v: &[(i64, i64)]) -> Ring {
        v.iter().map(|&(x, y)| p(x, y)).collect()
    }

    fn rect(x0: i64, y0: i64, x1: i64, y1: i64) -> Ring {
        ring(&[(x0, y0), (x1, y0), (x1, y1), (x0, y1)])
    }

    // The implementation before the bucket sorts and the narrow exact in-circle path, kept
    // to check that the current one gives identical triangulations.
    use crate::query::ring_area2;

    fn build_ref(polys: &[Polygon]) -> Result<Mesh> {
        let mut total = 0usize;
        for poly in polys {
            for ring in poly.rings() {
                for &p in ring.iter() {
                    if !p.in_range() {
                        return Err(Error::CoordinateOutOfRange(p));
                    }
                }
                total += ring.len();
            }
        }
        if total > u32::MAX as usize {
            return Err(Error::TooLarge);
        }
        // Directed edges, interior on the left.
        let mut raw: Vec<(Point, Point)> = Vec::with_capacity(total);
        let mut ring_buf: Vec<Point> = Vec::new();
        for poly in polys {
            for (ri, ring) in poly.rings().enumerate() {
                ring_buf.clear();
                for &p in ring.iter() {
                    if ring_buf.last() != Some(&p) {
                        ring_buf.push(p);
                    }
                }
                while ring_buf.len() > 1 && ring_buf.first() == ring_buf.last() {
                    ring_buf.pop();
                }
                if ring_buf.len() < 3 {
                    continue;
                }
                let a = ring_area2(&ring_buf);
                if a == 0 {
                    continue;
                }
                let flip = (ri == 0) != (a > 0);
                let n = ring_buf.len();
                for i in 0..n {
                    let (p, q) = (ring_buf[i], ring_buf[(i + 1) % n]);
                    raw.push(if flip { (q, p) } else { (p, q) });
                }
            }
        }
        // Split edges at vertices lying on their interiors (exactly), so that partially
        // shared edges become identical pieces that cancel below. Proper crossings mean
        // invalid input.
        let raw: Vec<(Point, Point)> = match crate::node::node_exact(&raw) {
            Ok(frags) => frags.into_iter().map(|f| (f.a, f.b)).collect(),
            Err(_) => return Err(INVALID),
        };
        let mut verts: Vec<Point> = Vec::with_capacity(raw.len());
        verts.extend(raw.iter().map(|e| e.0));
        verts.sort_unstable();
        verts.dedup();
        let id = |p: Point| verts.binary_search(&p).unwrap_or(0) as u32;
        let mut und: Vec<Edge> = raw
            .iter()
            .map(|&(a, b)| {
                let (ia, ib) = (id(a), id(b));
                if ia < ib {
                    Edge {
                        lo: ia,
                        hi: ib,
                        above: true,
                    }
                } else {
                    Edge {
                        lo: ib,
                        hi: ia,
                        above: false,
                    }
                }
            })
            .collect();
        drop(raw);
        und.sort_unstable_by_key(|e| (e.lo, e.hi));
        // Cancel opposite copies of the same segment; same-direction copies overlap.
        let mut edges: Vec<Edge> = Vec::with_capacity(und.len());
        let mut i = 0;
        while i < und.len() {
            let mut j = i;
            let mut net = 0i64;
            while j < und.len() && und[j].lo == und[i].lo && und[j].hi == und[i].hi {
                net += if und[j].above { 1 } else { -1 };
                j += 1;
            }
            match net {
                0 => {}
                1 | -1 => edges.push(Edge {
                    above: net > 0,
                    ..und[i]
                }),
                _ => return Err(INVALID),
            }
            i = j;
        }
        drop(und);
        let dir = |e: &Edge| sub(verts[e.hi as usize], verts[e.lo as usize]);
        edges.sort_by(|a, b| {
            a.lo.cmp(&b.lo)
                .then_with(|| cmp_dir_halfplane(dir(a), dir(b)))
        });
        let mut area2: i128 = 0;
        for e in &edges {
            let c = cross(verts[e.lo as usize], verts[e.hi as usize]);
            area2 += if e.above { c } else { -c };
        }
        Ok(Mesh {
            verts,
            edges,
            area2,
        })
    }

    fn finish_ref(verts: Vec<Point>, mut tris: Vec<[u32; 3]>) -> Triangulation {
        let mut used = vec![false; verts.len()];
        for t in &tris {
            for &i in t {
                used[i as usize] = true;
            }
        }
        let (vertices, remap) = if used.iter().all(|&u| u) {
            (verts, None)
        } else {
            let mut remap = vec![0u32; verts.len()];
            let mut out = Vec::with_capacity(verts.len());
            for (i, p) in verts.into_iter().enumerate() {
                if used[i] {
                    remap[i] = out.len() as u32;
                    out.push(p);
                }
            }
            (out, Some(remap))
        };
        for t in &mut tris {
            if let Some(r) = &remap {
                for i in t.iter_mut() {
                    *i = r[*i as usize];
                }
            }
            let k = (0..3).min_by_key(|&k| t[k]).unwrap_or(0);
            t.rotate_left(k);
        }
        tris.sort_unstable();
        Triangulation {
            vertices,
            triangles: tris,
        }
    }

    fn delaunay_flip_ref(pts: &[Point], tris: &mut [[u32; 3]]) {
        let nt = tris.len();
        // nbr[t][i]: triangle across the edge opposite corner i, NONE on the boundary.
        let mut nbr = vec![[NONE; 3]; nt];
        let mut half: Vec<(u32, u32, u32, u8)> = Vec::with_capacity(nt * 3);
        for (t, tri) in tris.iter().enumerate() {
            for i in 0..3 {
                let (a, b) = (tri[(i + 1) % 3], tri[(i + 2) % 3]);
                half.push((a.min(b), a.max(b), t as u32, i as u8));
            }
        }
        half.sort_unstable();
        let mut stack: Vec<(u32, u8)> = Vec::new();
        let mut i = 0;
        while i < half.len() {
            let mut j = i;
            while j < half.len() && half[j].0 == half[i].0 && half[j].1 == half[i].1 {
                j += 1;
            }
            if j - i == 2 {
                let (_, _, t, a) = half[i];
                let (_, _, u, b) = half[i + 1];
                nbr[t as usize][a as usize] = u;
                nbr[u as usize][b as usize] = t;
                stack.push((t, a));
            }
            i = j;
        }
        drop(half);
        stack.reverse();
        while let Some((t, i)) = stack.pop() {
            let (t, i) = (t as usize, i as usize);
            let u = nbr[t][i];
            if u == NONE {
                continue;
            }
            let u = u as usize;
            let Some(j) = (0..3).find(|&j| nbr[u][j] == t as u32) else {
                continue;
            };
            let (a, b, c) = (tris[t][i], tris[t][(i + 1) % 3], tris[t][(i + 2) % 3]);
            let d = tris[u][j];
            if tris[u][(j + 1) % 3] != c || tris[u][(j + 2) % 3] != b {
                continue;
            }
            let (pa, pb, pc, pd) = (
                pts[a as usize],
                pts[b as usize],
                pts[c as usize],
                pts[d as usize],
            );
            if incircle_ref(pa, pb, pc, pd) <= 0
                || orient(pa, pb, pd) <= 0
                || orient(pa, pd, pc) <= 0
            {
                continue;
            }
            let n_ab = nbr[t][(i + 2) % 3];
            let n_ca = nbr[t][(i + 1) % 3];
            let n_bd = nbr[u][(j + 1) % 3];
            let n_dc = nbr[u][(j + 2) % 3];
            tris[t] = [a, b, d];
            nbr[t] = [n_bd, u as u32, n_ab];
            tris[u] = [a, d, c];
            nbr[u] = [n_dc, n_ca, t as u32];
            if n_bd != NONE {
                for x in nbr[n_bd as usize].iter_mut() {
                    if *x == u as u32 {
                        *x = t as u32;
                    }
                }
            }
            if n_ca != NONE {
                for x in nbr[n_ca as usize].iter_mut() {
                    if *x == t as u32 {
                        *x = u as u32;
                    }
                }
            }
            stack.push((t as u32, 0));
            stack.push((t as u32, 2));
            stack.push((u as u32, 0));
            stack.push((u as u32, 1));
        }
    }

    /// Exact in-circle test: positive when `d` lies strictly inside the circumcircle of the
    /// counter-clockwise triangle `abc`, zero when co-circular, negative outside.
    fn incircle_ref(a: Point, b: Point, c: Point, d: Point) -> i32 {
        let adx = (a.x - d.x) as i128;
        let ady = (a.y - d.y) as i128;
        let bdx = (b.x - d.x) as i128;
        let bdy = (b.y - d.y) as i128;
        let cdx = (c.x - d.x) as i128;
        let cdy = (c.y - d.y) as i128;
        // Floating-point filter (Shewchuk's bound for the translated determinant). The
        // differences (|v| <= 2^41) are exact in f64.
        {
            let (fadx, fady, fbdx, fbdy, fcdx, fcdy) = (
                adx as f64, ady as f64, bdx as f64, bdy as f64, cdx as f64, cdy as f64,
            );
            let alift = fadx * fadx + fady * fady;
            let blift = fbdx * fbdx + fbdy * fbdy;
            let clift = fcdx * fcdx + fcdy * fcdy;
            let bc = fbdx * fcdy - fbdy * fcdx;
            let ca = fcdx * fady - fcdy * fadx;
            let ab = fadx * fbdy - fady * fbdx;
            let det = alift * bc + blift * ca + clift * ab;
            let perm = alift * ((fbdx * fcdy).abs() + (fbdy * fcdx).abs())
                + blift * ((fcdx * fady).abs() + (fcdy * fadx).abs())
                + clift * ((fadx * fbdy).abs() + (fady * fbdx).abs());
            let bound = 4e-15 * perm;
            if det > bound {
                return 1;
            }
            if det < -bound {
                return -1;
            }
        }
        let alift = adx * adx + ady * ady;
        let blift = bdx * bdx + bdy * bdy;
        let clift = cdx * cdx + cdy * cdy;
        let bc = bdx * cdy - bdy * cdx;
        let ca = cdx * ady - cdy * adx;
        let ab = adx * bdy - ady * bdx;
        let s = I256::mul(alift, bc)
            .add(I256::mul(blift, ca))
            .add(I256::mul(clift, ab));
        s.signum()
    }

    fn triangulate_ref(polys: &[Polygon], delaunay: bool) -> Result<Triangulation> {
        let mesh = build_ref(polys)?;
        let mut tris = mesh.sweep(CHUNK)?;
        if delaunay {
            delaunay_flip_ref(&mesh.verts, &mut tris);
        }
        Ok(finish_ref(mesh.verts, tris))
    }

    #[test]
    fn same_as_reference() {
        let mut s: u64 = 0xdead;
        let mut rnd = |m: i64| {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((s >> 33) as i64).rem_euclid(m)
        };
        let mut nonempty = 0;
        for it in 0..3000 {
            // A frame minus random rectangles and triangles (valid, many cocircular
            // points), at small and large scales; or raw rings (often invalid).
            let scale = [1i64, 1000, 1 << 32][it % 3];
            let frame = rect(0, 0, 40 * scale, 40 * scale);
            let k = 1 + rnd(12) as usize;
            let holes: Vec<Ring> = (0..k)
                .map(|_| {
                    let (x, y) = (rnd(36) * scale, rnd(36) * scale);
                    if rnd(3) == 0 {
                        ring(&[
                            (x, y),
                            (x + rnd(6) * scale + scale, y),
                            (x, y + rnd(6) * scale + scale),
                        ])
                    } else {
                        rect(x, y, x + (1 + rnd(5)) * scale, y + (1 + rnd(5)) * scale)
                    }
                })
                .collect();
            let polys: Vec<Polygon> = if it % 5 == 4 {
                (0..1 + rnd(3))
                    .map(|_| {
                        let v: Vec<(i64, i64)> = (0..3 + rnd(6))
                            .map(|_| (rnd(20) * scale, rnd(20) * scale))
                            .collect();
                        Polygon::new(ring(&v), vec![])
                    })
                    .collect()
            } else {
                crate::boolean::boolean(
                    crate::boolean::Op::Difference,
                    &frame,
                    &holes,
                    crate::boolean::FillRule::NonZero,
                )
                .unwrap()
            };
            for delaunay in [false, true] {
                let (a, b) = if delaunay && polys.len() == 1 {
                    (
                        triangulate_delaunay(&polys[0]),
                        triangulate_ref(&polys, true),
                    )
                } else if delaunay {
                    continue;
                } else {
                    (triangulate_set(&polys), triangulate_ref(&polys, false))
                };
                nonempty += a.as_ref().is_ok_and(|t| !t.triangles.is_empty()) as usize;
                assert_eq!(a, b, "{polys:?}");
            }
        }
        assert!(nonempty > 1000);
    }

    /// Canonically oriented, duplicate-free directed boundary edges of a polygon set.
    fn boundary(polys: &[Polygon]) -> Vec<(Point, Point)> {
        let mut out = Vec::new();
        for poly in polys {
            for (ri, r) in poly.rings().enumerate() {
                let mut v: Vec<Point> = Vec::new();
                for &q in r.iter() {
                    if v.last() != Some(&q) {
                        v.push(q);
                    }
                }
                while v.len() > 1 && v.first() == v.last() {
                    v.pop();
                }
                let a = ring_area2(&v);
                if v.len() < 3 || a == 0 {
                    continue;
                }
                if (ri == 0) != (a > 0) {
                    v.reverse();
                }
                for i in 0..v.len() {
                    out.push((v[i], v[(i + 1) % v.len()]));
                }
            }
        }
        out
    }

    /// Checks every documented invariant of a triangulation of a valid polygon set.
    fn check(polys: &[Polygon], t: &Triangulation) {
        assert!(
            t.vertices.windows(2).all(|w| w[0] < w[1]),
            "vertices sorted"
        );
        assert!(
            t.triangles.windows(2).all(|w| w[0] < w[1]),
            "triangles sorted"
        );
        let mut half: HashMap<(u32, u32), u32> = HashMap::new();
        for (i, tri) in t.triangles.iter().enumerate() {
            assert!(
                tri[0] < tri[1] && tri[0] < tri[2],
                "canonical start {tri:?}"
            );
            let [a, b, c] = t.triangle(i);
            assert!(orient(a, b, c) > 0, "bad triangle {a:?} {b:?} {c:?}");
            for k in 0..3 {
                *half.entry((tri[k], tri[(k + 1) % 3])).or_default() += 1;
            }
        }
        let bnd_edges = boundary(polys);
        let expect: i128 = bnd_edges.iter().map(|&(a, b)| cross(a, b)).sum();
        assert_eq!(t.area2(), expect, "area");
        let id = |q: Point| t.vertices.binary_search(&q).expect("vertex present") as u32;
        let mut bnd: HashMap<(u32, u32), u32> = HashMap::new();
        for &(a, b) in &bnd_edges {
            *bnd.entry((id(a), id(b))).or_default() += 1;
        }
        for (&(a, b), &n) in &bnd {
            assert_eq!(n, 1);
            assert_eq!(half.get(&(a, b)), Some(&1), "boundary edge missing");
            assert_eq!(half.get(&(b, a)), None, "triangle outside boundary edge");
        }
        for (&(a, b), &n) in &half {
            assert_eq!(n, 1, "edge used twice in the same direction");
            if !bnd.contains_key(&(a, b)) {
                assert_eq!(half.get(&(b, a)), Some(&1), "interior edge not shared");
            }
        }
        let mut all: Vec<Point> = bnd_edges.iter().map(|e| e.0).collect();
        all.sort_unstable();
        all.dedup();
        assert_eq!(all, t.vertices, "vertex set");
    }

    /// Every non-boundary edge is locally Delaunay.
    fn check_delaunay(t: &Triangulation) {
        let mut opp: HashMap<(u32, u32), u32> = HashMap::new();
        for tri in &t.triangles {
            for k in 0..3 {
                opp.insert((tri[k], tri[(k + 1) % 3]), tri[(k + 2) % 3]);
            }
        }
        let v = |i: u32| t.vertices[i as usize];
        for tri in &t.triangles {
            for k in 0..3 {
                let (a, b, c) = (tri[k], tri[(k + 1) % 3], tri[(k + 2) % 3]);
                if let Some(&d) = opp.get(&(b, a)) {
                    assert!(incircle(v(a), v(b), v(c), v(d)) <= 0, "not Delaunay");
                }
            }
        }
    }

    fn run(poly: &Polygon) -> Triangulation {
        let t = triangulate(poly).expect("triangulate");
        check(core::slice::from_ref(poly), &t);
        // Tiny status chunks (exercising chunk splits and merges) give the same result.
        let mesh = Mesh::build(core::slice::from_ref(poly)).unwrap();
        let tiny = finish(mesh.verts.clone(), mesh.sweep(1).unwrap());
        assert_eq!(tiny, t);
        let d = triangulate_delaunay(poly).expect("delaunay");
        check(core::slice::from_ref(poly), &d);
        check_delaunay(&d);
        // Deterministic and independent of ring start, orientation and hole order.
        assert_eq!(triangulate(poly).unwrap(), t);
        let mut q = poly.clone();
        q.outer.rotate_left(1);
        q.outer.reverse();
        for h in &mut q.holes {
            let k = 2.min(h.len());
            h.rotate_left(k);
        }
        q.holes.reverse();
        assert_eq!(triangulate(&q).unwrap(), t);
        t
    }

    #[test]
    fn square_and_triangle() {
        let t = run(&rect(0, 0, 10, 10).into());
        assert_eq!(t.triangles.len(), 2);
        let t = run(&ring(&[(0, 0), (10, 0), (3, 7)]).into());
        assert_eq!(t.triangles, vec![[0, 2, 1]]);
    }

    #[test]
    fn concave() {
        run(&ring(&[(0, 0), (10, 0), (10, 2), (2, 2), (2, 10), (0, 10)]).into());
        run(&ring(&[(0, 0), (10, 5), (0, 10), (3, 5)]).into());
        run(&ring(&[(0, 0), (5, 3), (10, 0), (10, 10), (5, 7), (0, 10)]).into());
        let mut v = Vec::new();
        for i in 0..16 {
            let a = i as f64 * core::f64::consts::PI / 8.0;
            let r = if i % 2 == 0 { 1000.0 } else { 300.0 };
            v.push(((a.cos() * r) as i64, (a.sin() * r) as i64));
        }
        run(&ring(&v).into());
    }

    #[test]
    fn combs() {
        // Teeth pointing up.
        let mut up = vec![(0, 0), (100, 0)];
        for i in (0..10).rev() {
            up.extend([(i * 10 + 10, 50), (i * 10 + 5, 50), (i * 10 + 5, 10)]);
            if i > 0 {
                up.push((i * 10, 10));
            }
        }
        up.push((0, 50));
        run(&ring(&up).into());
        // Slots cut in from the right (merge vertices), mirrored (split vertices), rotated.
        let mut right = vec![(0, 0), (100, 0)];
        for i in 0..9 {
            let y = i * 10;
            right.extend([(100, y + 5), (20, y + 5), (20, y + 10), (100, y + 10)]);
        }
        right.extend([(100, 100), (0, 100)]);
        run(&ring(&right).into());
        let left: Vec<(i64, i64)> = right.iter().map(|&(x, y)| (100 - x, y)).collect();
        run(&ring(&left).into());
        let rot: Vec<(i64, i64)> = right.iter().map(|&(x, y)| (y, x)).collect();
        run(&ring(&rot).into());
        let skew: Vec<(i64, i64)> = right.iter().map(|&(x, y)| (x + y / 3, y - x / 7)).collect();
        run(&ring(&skew).into());
    }

    #[test]
    fn spiral() {
        // A thick square spiral built as a union of axis-aligned bars.
        let dirs = [(1, 0), (0, 1), (-1, 0), (0, -1)];
        let (mut x, mut y) = (0i64, 0i64);
        let mut bars = Vec::new();
        for k in 0..24 {
            let len = 6 * (k as i64 / 2 + 1);
            let (dx, dy) = dirs[k % 4];
            let (nx, ny) = (x + dx * len, y + dy * len);
            bars.push(rect(
                x.min(nx) - 1,
                y.min(ny) - 1,
                x.max(nx) + 1,
                y.max(ny) + 1,
            ));
            (x, y) = (nx, ny);
        }
        let polys = crate::union_all(&bars, crate::FillRule::NonZero).unwrap();
        assert_eq!(polys.len(), 1);
        run(&polys[0]);
        // Its complement in a box: a spiral with a hole.
        let boxed = crate::boolean(
            crate::Op::Difference,
            &rect(-100, -100, 100, 100),
            &bars,
            crate::FillRule::NonZero,
        )
        .unwrap();
        for poly in &boxed {
            run(poly);
        }
    }

    #[test]
    fn holes() {
        let poly = Polygon::new(
            rect(0, 0, 100, 100),
            vec![
                rect(10, 10, 20, 20),
                rect(30, 10, 40, 90),
                rect(50, 50, 90, 60),
            ],
        );
        run(&poly);
        // Diamond holes: split vertices strictly inside regions.
        let mut hs = Vec::new();
        for i in 0..5 {
            for j in 0..5 {
                let (cx, cy) = (10 + 20 * i, 10 + 20 * j + i);
                hs.push(ring(&[
                    (cx - 5, cy),
                    (cx, cy - 5),
                    (cx + 5, cy),
                    (cx, cy + 5),
                ]));
            }
        }
        run(&Polygon::new(rect(0, 0, 110, 110), hs));
    }

    #[test]
    fn touching_rings() {
        // Hole touching the outer ring at a vertex.
        let outer = ring(&[(0, 0), (10, 0), (10, 10), (0, 10), (0, 5)]);
        let hole = ring(&[(0, 5), (3, 7), (5, 5), (3, 3)]);
        run(&Polygon::new(outer.clone(), vec![hole.clone()]));
        // Two holes sharing a vertex, one of them touching the outer ring.
        let h2 = ring(&[(5, 5), (7, 8), (9, 5), (7, 2)]);
        run(&Polygon::new(outer, vec![hole, h2]));
        // Two triangular holes meeting at a point (both orientations of the pair).
        let c = (50, 50);
        let quads = [
            ring(&[c, (60, 40), (60, 60)]),
            ring(&[c, (60, 60), (40, 60)]),
            ring(&[c, (40, 60), (40, 40)]),
            ring(&[c, (40, 40), (60, 40)]),
        ];
        for k in 0..2 {
            let hs = vec![quads[k].clone(), quads[k + 2].clone()];
            run(&Polygon::new(rect(0, 0, 100, 100), hs));
        }
        // Hole touching the outer ring at the hole's leftmost vertex.
        let outer = ring(&[(0, 0), (20, 0), (20, 20), (0, 20), (0, 10)]);
        let hole = ring(&[(0, 10), (10, 15), (12, 10), (10, 5)]);
        run(&Polygon::new(outer, vec![hole]));
        // Hole touching the outer ring at a reflex outer vertex (from above and below).
        let outer = ring(&[(0, 0), (20, 0), (20, 20), (10, 12), (0, 20)]);
        run(&Polygon::new(
            outer,
            vec![ring(&[(10, 12), (8, 6), (12, 6)])],
        ));
        let outer = ring(&[(0, 0), (10, 8), (20, 0), (20, 20), (0, 20)]);
        run(&Polygon::new(
            outer,
            vec![ring(&[(10, 8), (8, 14), (12, 14)])],
        ));
        // Hole touching the outer ring at the hole's rightmost vertex.
        let outer = ring(&[(0, 0), (20, 0), (20, 10), (20, 20), (0, 20)]);
        run(&Polygon::new(
            outer,
            vec![ring(&[(20, 10), (10, 5), (8, 10), (10, 15)])],
        ));
        // Chain of holes touching each other.
        let mut hs = Vec::new();
        for i in 0..8 {
            let x = 10 + 10 * i;
            hs.push(ring(&[(x, 50), (x + 5, 45), (x + 10, 50), (x + 5, 55)]));
        }
        run(&Polygon::new(rect(0, 0, 100, 100), hs));
        // Vertical chain.
        let mut hs = Vec::new();
        for i in 0..8 {
            let y = 10 + 10 * i;
            hs.push(ring(&[(50, y), (55, y + 5), (50, y + 10), (45, y + 5)]));
        }
        run(&Polygon::new(rect(0, 0, 100, 100), hs));
    }

    #[test]
    fn polygon_sets() {
        // Squares touching at corners, an island inside a hole.
        let polys = vec![
            Polygon::from(rect(0, 0, 10, 10)),
            Polygon::from(rect(10, 10, 18, 18)),
            Polygon::new(rect(20, 0, 60, 40), vec![rect(25, 5, 55, 35)]),
            Polygon::from(rect(30, 10, 50, 30)),
        ];
        let t = triangulate_set(&polys).unwrap();
        check(&polys, &t);
        assert_eq!(t.vertices.iter().filter(|&&q| q == p(10, 10)).count(), 1);
        assert!(triangulate_set(&[]).unwrap().triangles.is_empty());
    }

    #[test]
    fn collinear_vertices() {
        let r = ring(&[
            (0, 0),
            (5, 0),
            (10, 0),
            (10, 5),
            (10, 10),
            (5, 10),
            (0, 10),
            (0, 5),
        ]);
        let t = run(&r.clone().into());
        assert_eq!(t.vertices.len(), 8);
        assert_eq!(t.triangles.len(), 6);
        let h = ring(&[
            (2, 2),
            (2, 4),
            (2, 6),
            (4, 6),
            (6, 6),
            (6, 4),
            (6, 2),
            (4, 2),
        ]);
        run(&Polygon::new(r, vec![h]));
        let mut v: Vec<(i64, i64)> = (0..20).map(|i| (i, 0)).collect();
        v.extend((0..20).map(|i| (20 - i, i)));
        run(&ring(&v).into());
        let mut v: Vec<(i64, i64)> = (0..20).map(|i| (0, -i)).collect();
        v.extend((0..20).map(|i| (i, -20)));
        v.push((20, 0));
        run(&ring(&v).into());
    }

    #[test]
    fn t_junction_split() {
        // Hole vertex in the interior of an outer edge.
        let poly = Polygon::new(rect(0, 0, 10, 10), vec![ring(&[(0, 5), (4, 7), (4, 3)])]);
        let t = triangulate(&poly).unwrap();
        assert_eq!(t.area2(), 200 - 16);
        assert!(t.vertices.contains(&p(0, 5)));
        for i in 0..t.triangles.len() {
            let [a, b, c] = t.triangle(i);
            assert!(orient(a, b, c) > 0);
        }
    }

    #[test]
    fn degenerate_inputs() {
        let tri = |r: Ring| triangulate(&r.into());
        assert!(tri(Ring::new()).unwrap().triangles.is_empty());
        assert!(tri(ring(&[(0, 0), (1, 1)])).unwrap().triangles.is_empty());
        assert!(
            tri(ring(&[(0, 0), (5, 5), (10, 10)]))
                .unwrap()
                .triangles
                .is_empty()
        );
        let dup = ring(&[(0, 0), (0, 0), (10, 0), (10, 0), (0, 10), (0, 0)]);
        assert_eq!(tri(dup).unwrap().triangles.len(), 1);
        // Spike (back-and-forth edge) is cancelled; its tip is dropped.
        let t = tri(ring(&[(0, 0), (10, 0), (20, 0), (10, 0), (10, 10)])).unwrap();
        assert_eq!(t.area2(), 100);
        assert!(!t.vertices.contains(&p(20, 0)));
        let big = crate::MAX_COORD + 1;
        assert_eq!(
            tri(ring(&[(0, 0), (big, 0), (0, 1)])),
            Err(Error::CoordinateOutOfRange(p(big, 0)))
        );
        assert_eq!(
            tri(ring(&[(0, 0), (1, 0), (0, -big)])),
            Err(Error::CoordinateOutOfRange(p(0, -big)))
        );
        // Invalid input: errors, never panics.
        assert!(tri(ring(&[(0, 0), (10, 10), (10, 0), (0, 20)])).is_err());
        let bad = [
            Polygon::new(rect(0, 0, 10, 10), vec![rect(20, 0, 30, 10)]),
            Polygon::new(rect(0, 0, 10, 10), vec![rect(2, 2, 6, 6), rect(4, 4, 8, 8)]),
            Polygon::new(rect(0, 0, 10, 10), vec![rect(5, 2, 15, 8)]),
            Polygon::new(rect(0, 0, 10, 10), vec![rect(2, 2, 6, 6), rect(2, 2, 6, 6)]),
        ];
        for poly in &bad {
            assert!(triangulate(poly).is_err(), "{poly:?}");
            assert!(triangulate_delaunay(poly).is_err(), "{poly:?}");
        }
        // A hole equal to the outer ring cancels it: empty result.
        let gone = Polygon::new(rect(0, 0, 10, 10), vec![rect(0, 0, 10, 10)]);
        assert_eq!(triangulate(&gone).unwrap(), Triangulation::default());
        let polys = [
            Polygon::from(rect(0, 0, 10, 10)),
            Polygon::from(rect(5, 5, 15, 15)),
        ];
        assert!(triangulate_set(&polys).is_err());
        // Zero-area hole is ignored.
        let zh = Polygon::new(rect(0, 0, 10, 10), vec![ring(&[(2, 2), (4, 4), (6, 6)])]);
        assert_eq!(triangulate(&zh).unwrap().area2(), 200);
        // Polygons sharing an edge: the shared edge cancels.
        let polys = [
            Polygon::from(rect(0, 0, 10, 10)),
            Polygon::from(rect(10, 0, 20, 10)),
        ];
        assert_eq!(triangulate_set(&polys).unwrap().area2(), 400);
    }

    #[test]
    fn extreme_coordinates() {
        let m = crate::MAX_COORD;
        let poly = Polygon::new(
            rect(-m, -m, m, m),
            vec![ring(&[(-m + 1, 0), (0, m - 1), (m - 1, 0), (0, -m + 1)])],
        );
        run(&poly);
        // Nearly co-circular points at full range exercise the exact in-circle path.
        let r = ring(&[(-m, -m), (m, -m), (m, m), (-m, m), (-m + 1, 0)]);
        run(&r.into());
    }

    #[test]
    fn incircle_exact() {
        let m = crate::MAX_COORD;
        assert_eq!(incircle(p(-m, -m), p(m, -m), p(m, m), p(-m, m)), 0);
        assert_eq!(incircle(p(-m, -m), p(m, -m), p(m, m), p(-m + 1, m - 1)), 1);
        assert_eq!(incircle(p(-m, -m), p(m, -m), p(m - 1, m), p(-m, m)), -1);
        // Agreement with a direct i128 evaluation on small coordinates.
        let mut s = 12345u64;
        let mut rnd = || {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((s >> 33) % 41) as i64 - 20
        };
        for _ in 0..20000 {
            let q: Vec<Point> = (0..4).map(|_| p(rnd(), rnd())).collect();
            let (a, b, c, d) = (q[0], q[1], q[2], q[3]);
            let f = |u: Point| ((u.x - d.x) as i128, (u.y - d.y) as i128);
            let ((ax, ay), (bx, by), (cx, cy)) = (f(a), f(b), f(c));
            let det = (ax * ax + ay * ay) * (bx * cy - by * cx)
                + (bx * bx + by * by) * (cx * ay - cy * ax)
                + (cx * cx + cy * cy) * (ax * by - ay * bx);
            assert_eq!(incircle(a, b, c, d), det.signum() as i32);
        }
        let big = (1i128 << 100) + 12345;
        let r = I256::mul(big, -big);
        assert_eq!(r.signum(), -1);
        assert_eq!(r.add(I256::mul(big, big)), I256 { hi: 0, lo: 0 });
        assert_eq!(I256::mul(-3, 5).add(I256::mul(3, 5)).signum(), 0);
        assert_eq!(I256::mul(1 << 90, 1 << 90), I256 { hi: 1 << 52, lo: 0 });
    }

    proptest::proptest! {
        /// The chunked status gives the same result for any chunk length.
        #[test]
        fn chunk_length_irrelevant(
            rs in proptest::collection::vec(
                proptest::collection::vec((-30i64..=30, -30i64..=30), 0..12),
                0..8,
            )
        ) {
            let rings: Vec<Ring> = rs.iter().map(|r| ring(r)).collect();
            let polys = crate::union_all(&rings, crate::FillRule::EvenOdd).unwrap();
            let mesh = Mesh::build(&polys).unwrap();
            let a = mesh.sweep(CHUNK).unwrap();
            for chunk in [1, 2, 3] {
                proptest::prop_assert_eq!(&mesh.sweep(chunk).unwrap(), &a);
            }
        }
    }

    /// A zone fill: a large outline minus 2000 overlapping 50-gon obstacles.
    fn zone() -> Polygon {
        let mut s = 7u64;
        let mut rnd = |n: u64| {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (s >> 33) % n
        };
        let mut obstacles: Vec<Ring> = Vec::new();
        for i in 0..45 {
            for j in 0..45 {
                if obstacles.len() == 2000 {
                    break;
                }
                let cx = 1000 + i * 2000 + rnd(600) as i64;
                let cy = 1000 + j * 2000 + rnd(600) as i64;
                let r = 300.0 + rnd(400) as f64;
                obstacles.push(
                    (0..50)
                        .map(|t| {
                            let a = t as f64 * core::f64::consts::TAU / 50.0;
                            p(cx + (a.cos() * r) as i64, cy + (a.sin() * r) as i64)
                        })
                        .collect(),
                );
            }
        }
        let outline = rect(0, 0, 91000, 91000);
        let mut res = crate::boolean(
            crate::Op::Difference,
            &outline,
            &obstacles,
            crate::FillRule::NonZero,
        )
        .unwrap();
        res.sort_by_key(|q| core::cmp::Reverse(q.vertex_count()));
        res.swap_remove(0)
    }

    #[test]
    fn zone_fill_performance() {
        let poly = zone();
        assert!(poly.holes.len() > 1000, "{} holes", poly.holes.len());
        let n = poly.vertex_count();
        let t0 = std::time::Instant::now();
        let t = triangulate(&poly).unwrap();
        let dt = t0.elapsed();
        let t1 = std::time::Instant::now();
        let d = triangulate_delaunay(&poly).unwrap();
        let dd = t1.elapsed();
        eprintln!(
            "zone: {} holes, {n} vertices, {} triangles: sweep {dt:?}, delaunay {dd:?}",
            poly.holes.len(),
            t.triangles.len()
        );
        check(core::slice::from_ref(&poly), &t);
        check(core::slice::from_ref(&poly), &d);
        let mesh = Mesh::build(core::slice::from_ref(&poly)).unwrap();
        assert_eq!(finish(mesh.verts.clone(), mesh.sweep(2).unwrap()), t);
        check_delaunay(&d);
        if !cfg!(debug_assertions) {
            assert!(dt.as_secs_f64() < 0.5, "sweep too slow: {dt:?}");
            assert!(dd.as_secs_f64() < 1.0, "delaunay too slow: {dd:?}");
        }
    }

    #[test]
    fn wide_status_performance() {
        // 50 000 thin slots stacked vertically, starting at scattered x: the sweep status
        // grows to 100 000 edges and every insertion lands at a random position.
        let k = 50_000i64;
        let holes: Vec<Ring> = (0..k)
            .map(|i| {
                let (x, y) = (10 + (i * 7919) % 60_000, 10 + 4 * i);
                rect(x, y, 100_000, y + 2)
            })
            .collect();
        let poly = Polygon::new(rect(0, 0, 200_000, 20 + 4 * k), holes);
        let t0 = std::time::Instant::now();
        let t = triangulate(&poly).unwrap();
        let dt = t0.elapsed();
        eprintln!("wide: {} triangles in {dt:?}", t.triangles.len());
        check(core::slice::from_ref(&poly), &t);
        let mesh = Mesh::build(core::slice::from_ref(&poly)).unwrap();
        assert_eq!(finish(mesh.verts.clone(), mesh.sweep(3).unwrap()), t);
        if !cfg!(debug_assertions) {
            assert!(dt.as_secs_f64() < 0.5, "sweep too slow: {dt:?}");
        }
    }
}
