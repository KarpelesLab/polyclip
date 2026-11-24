//! Results of independent clusters of rings, kept from one operation to the next.
//!
//! A boolean splits its input into clusters of rings that cannot interact (see
//! `cluster.rs`) and computes each on its own. What a cluster gives depends only on its
//! edges, on the winding of the other rings around it and on the fill rules and operation,
//! and it moves with the cluster: every predicate is exact and looks at differences of
//! coordinates only, and new vertices are exact roundings of crossings. So the result of a
//! cluster can be reused for any cluster with the same edges up to a translation, within one
//! operation (done always) and, with a [`Memo`], across operations: a zone refilled after a
//! track moved has most of its clusters as they were.

use crate::arrangement::InEdge;
use crate::assemble::RawRing;
use crate::cluster::Part;
use crate::geom::Point;
use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex, MutexGuard};

/// One step of the signature hash.
#[inline]
pub(crate) fn mix(h: u64, w: u64) -> u64 {
    (h.rotate_left(5) ^ w).wrapping_mul(0x517c_c1b7_2722_0a95)
}

/// Hash of the rings of a cluster relative to `origin`, with the winding `base` of the
/// other clusters around it.
pub(crate) fn signature<'a>(
    rings: impl Iterator<Item = &'a [InEdge]>,
    origin: Point,
    base: [i32; 2],
) -> u64 {
    let mut h = mix(mix(0x9e37_79b9_7f4a_7c15, base[0] as u64), base[1] as u64);
    for r in rings {
        h = mix(h, r.len() as u64);
        for e in r {
            h = mix(h, e.a.x.wrapping_sub(origin.x) as u64);
            h = mix(h, e.a.y.wrapping_sub(origin.y) as u64);
            h = mix(h, e.tag ^ ((e.operand as u64) << 56));
        }
    }
    h
}

/// Whether two clusters have the same edges up to the translation taking `oa` to `ob`.
pub(crate) fn same_rings<'a>(
    a: impl Iterator<Item = &'a [InEdge]>,
    oa: Point,
    mut b: impl Iterator<Item = &'a [InEdge]>,
    ob: Point,
) -> bool {
    let (dx, dy) = (ob.x.wrapping_sub(oa.x), ob.y.wrapping_sub(oa.y));
    let moved = |p: Point, q: Point| p.x.wrapping_add(dx) == q.x && p.y.wrapping_add(dy) == q.y;
    for ra in a {
        let Some(rb) = b.next() else {
            return false;
        };
        if ra.len() != rb.len()
            || !ra.iter().zip(rb).all(|(e, f)| {
                e.tag == f.tag && e.operand == f.operand && moved(e.a, f.a) && moved(e.b, f.b)
            })
        {
            return false;
        }
    }
    b.next().is_none()
}

/// Points relative to the first vertex of a cluster, in 32 bits each way where all fit.
enum Pts {
    Narrow(Vec<[i32; 2]>),
    Wide(Vec<Point>),
}

impl Pts {
    fn new(pts: Vec<Point>) -> Pts {
        let narrow = |p: &Point| i32::try_from(p.x).is_ok() && i32::try_from(p.y).is_ok();
        if pts.iter().all(narrow) {
            Pts::Narrow(pts.iter().map(|p| [p.x as i32, p.y as i32]).collect())
        } else {
            Pts::Wide(pts)
        }
    }

    #[inline]
    fn get(&self, i: usize) -> Point {
        match self {
            Pts::Narrow(v) => Point::new(v[i][0] as i64, v[i][1] as i64),
            Pts::Wide(v) => v[i],
        }
    }

    fn bytes(&self) -> usize {
        match self {
            Pts::Narrow(v) => 8 * v.len(),
            Pts::Wide(v) => 16 * v.len(),
        }
    }
}

/// Tags of edges; nothing when all are zero.
fn pack_tags(tags: Vec<u64>) -> Vec<u64> {
    if tags.iter().all(|&t| t == 0) {
        Vec::new()
    } else {
        tags
    }
}

/// A cluster and its result, both relative to the cluster's first vertex.
struct Entry {
    /// Fill rules, operation and collinear flag of the boolean.
    code: u32,
    base: [i32; 2],
    /// Length and operand of every ring.
    rings: Vec<(u32, u8)>,
    /// Start of every edge (rings are closed chains).
    pts: Pts,
    /// Tag of every edge.
    tags: Vec<u64>,
    /// The result: length of every ring, their points and tags one ring after another.
    out_rings: Vec<u32>,
    out_pts: Pts,
    out_tags: Vec<u64>,
    is_hole: Vec<bool>,
    parent: Vec<Option<u32>>,
    bytes: usize,
}

impl Entry {
    /// `None` when the rings are not closed chains of one operand each.
    fn new<'a>(
        code: u32,
        base: [i32; 2],
        cluster: impl Iterator<Item = &'a [InEdge]>,
        origin: Point,
        part: &Part,
    ) -> Option<Entry> {
        let rel = |p: Point| Point::new(p.x.wrapping_sub(origin.x), p.y.wrapping_sub(origin.y));
        let mut rings = Vec::new();
        let mut pts = Vec::new();
        let mut tags = Vec::new();
        for r in cluster {
            let n = r.len();
            for (i, e) in r.iter().enumerate() {
                if e.operand != r[0].operand || e.b != r[if i + 1 == n { 0 } else { i + 1 }].a {
                    return None;
                }
                pts.push(rel(e.a));
                tags.push(e.tag);
            }
            rings.push((u32::try_from(n).ok()?, r.first()?.operand));
        }
        let (pts, tags) = (Pts::new(pts), pack_tags(tags));
        let out_rings: Vec<u32> = part.rings.iter().map(|r| r.pts.len() as u32).collect();
        let out_pts = Pts::new(
            part.rings
                .iter()
                .flat_map(|r| r.pts.iter().map(|&p| rel(p)))
                .collect(),
        );
        let out_tags = pack_tags(
            part.rings
                .iter()
                .flat_map(|r| r.tags.iter().copied())
                .collect(),
        );
        let bytes = 160
            + 8 * rings.len()
            + pts.bytes()
            + 8 * tags.len()
            + 16 * out_rings.len()
            + out_pts.bytes()
            + 8 * out_tags.len();
        Some(Entry {
            code,
            base,
            rings,
            pts,
            tags,
            out_rings,
            out_pts,
            out_tags,
            is_hole: part.is_hole.clone(),
            parent: part.parent.clone(),
            bytes,
        })
    }

    fn matches<'a>(
        &self,
        code: u32,
        base: [i32; 2],
        mut cluster: impl Iterator<Item = &'a [InEdge]>,
        origin: Point,
    ) -> bool {
        if self.code != code || self.base != base {
            return false;
        }
        let mut k = 0;
        for &(n, operand) in &self.rings {
            let Some(r) = cluster.next() else {
                return false;
            };
            if r.len() != n as usize {
                return false;
            }
            for (i, e) in r.iter().enumerate() {
                let p = self.pts.get(k);
                if e.operand != operand
                    || e.a.x.wrapping_sub(origin.x) != p.x
                    || e.a.y.wrapping_sub(origin.y) != p.y
                    || e.b != r[if i + 1 == r.len() { 0 } else { i + 1 }].a
                    || e.tag != self.tags.get(k).copied().unwrap_or(0)
                {
                    return false;
                }
                k += 1;
            }
        }
        cluster.next().is_none()
    }

    /// The result, with the cluster's first vertex at `origin`.
    fn part(&self, origin: Point) -> Part {
        let mut k = 0;
        let rings = self
            .out_rings
            .iter()
            .map(|&n| {
                let n = n as usize;
                let pts = (k..k + n)
                    .map(|i| {
                        let p = self.out_pts.get(i);
                        Point::new(p.x.wrapping_add(origin.x), p.y.wrapping_add(origin.y))
                    })
                    .collect();
                let tags = if self.out_tags.is_empty() {
                    vec![0; n]
                } else {
                    self.out_tags[k..k + n].to_vec()
                };
                k += n;
                RawRing { pts, tags }
            })
            .collect();
        Part {
            rings,
            is_hole: self.is_hole.clone(),
            parent: self.parent.clone(),
        }
    }
}

/// The hashes are well mixed already.
#[derive(Default)]
struct Identity(u64);

impl Hasher for Identity {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = (self.0 << 8) | b as u64;
        }
    }
    fn write_u64(&mut self, v: u64) {
        self.0 = v;
    }
}

struct Slot {
    entry: Arc<Entry>,
    used: u64,
}

#[derive(Default)]
struct Inner {
    map: HashMap<u64, Vec<Slot>, BuildHasherDefault<Identity>>,
    bytes: usize,
    tick: u64,
    hits: u64,
    misses: u64,
}

struct Shared {
    max_bytes: usize,
    inner: Mutex<Inner>,
}

/// Results of earlier booleans and offsets, reused where a later one has the same rings.
///
/// A boolean computes groups of rings that cannot interact (*clusters*) independently. A
/// `Memo` keeps the results of clusters, so that an operation on input that mostly repeats
/// an earlier one (a copper zone filled again after a few tracks moved, the same pad pattern
/// at another place) only computes the clusters that changed. A cluster is recognized by its
/// exact edges, up to a translation, together with the fill rules and the operation.
///
/// **The output of every operation is the same with and without a memo**; only the time
/// differs. Inputs too small to be split into clusters do not use it.
///
/// A `Memo` is a handle: clones share the same store, and it can be used from several
/// threads at once. It holds at most the number of bytes it was created with (approximately)
/// and drops the results used least recently beyond that.
///
/// ```
/// use polyclip::{ArcTol, Memo, Ring, Side, opening};
/// let memo = Memo::new(16 << 20);
/// let pads: Vec<Ring> = (0..200)
///     .map(|i| {
///         let (x, y) = (1000 * (i % 20), 1000 * (i / 20));
///         Ring::from([(x, y), (x + 400, y), (x + 400, y + 300), (x + 200, y + 400), (x, y + 300)])
///     })
///     .collect();
/// let tol = ArcTol::new(2, Side::Inside);
/// let first = memo.opening(&pads, 50, tol).unwrap();
/// let again = memo.opening(&pads, 50, tol).unwrap();
/// assert_eq!(first, again);
/// assert_eq!(first, opening(&pads, 50, tol).unwrap());
/// ```
#[derive(Clone)]
pub struct Memo(Arc<Shared>);

/// What a [`Memo`] holds and how often it was of use.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MemoStats {
    /// Clusters kept.
    pub entries: usize,
    /// Approximate size of what is kept, in bytes.
    pub bytes: usize,
    /// Clusters whose result was found.
    pub hits: u64,
    /// Clusters that were looked up, not found and computed.
    pub misses: u64,
}

impl Default for Memo {
    /// A memo of 64 MiB.
    fn default() -> Self {
        Memo::new(64 << 20)
    }
}

impl core::fmt::Debug for Memo {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("Memo")
            .field("max_bytes", &self.0.max_bytes)
            .field("stats", &self.stats())
            .finish()
    }
}

impl Memo {
    /// Creates an empty memo holding at most about `max_bytes` bytes of results.
    pub fn new(max_bytes: usize) -> Self {
        Memo(Arc::new(Shared {
            max_bytes,
            inner: Mutex::new(Inner::default()),
        }))
    }

    fn lock(&self) -> MutexGuard<'_, Inner> {
        // Whatever a panicking thread left is a consistent store of exact results.
        self.0.inner.lock().unwrap_or_else(|e| e.into_inner())
    }

    /// Forgets every result (and the counts of [`stats`](Self::stats)).
    pub fn clear(&self) {
        *self.lock() = Inner::default();
    }

    /// Number and size of the results kept, and how many lookups found one.
    pub fn stats(&self) -> MemoStats {
        let g = self.lock();
        MemoStats {
            entries: g.map.values().map(Vec::len).sum(),
            bytes: g.bytes,
            hits: g.hits,
            misses: g.misses,
        }
    }

    /// The result kept for the cluster with signature `sig`, moved to `origin`.
    pub(crate) fn get<'a>(
        &self,
        sig: u64,
        code: u32,
        base: [i32; 2],
        cluster: impl Iterator<Item = &'a [InEdge]> + Clone,
        origin: Point,
    ) -> Option<Part> {
        let key = mix(sig, code as u64);
        let found: Vec<Arc<Entry>> = {
            let mut g = self.lock();
            g.tick += 1;
            let tick = g.tick;
            match g.map.get_mut(&key) {
                Some(slots) => slots
                    .iter_mut()
                    .map(|s| {
                        s.used = tick;
                        s.entry.clone()
                    })
                    .collect(),
                None => Vec::new(),
            }
        };
        let hit = found
            .iter()
            .find(|e| e.matches(code, base, cluster.clone(), origin));
        let mut g = self.lock();
        match hit {
            Some(e) => {
                g.hits += 1;
                drop(g);
                Some(e.part(origin))
            }
            None => {
                g.misses += 1;
                None
            }
        }
    }

    /// Keeps `part`, the result of the cluster with signature `sig` at `origin`.
    pub(crate) fn put<'a>(
        &self,
        sig: u64,
        code: u32,
        base: [i32; 2],
        cluster: impl Iterator<Item = &'a [InEdge]>,
        origin: Point,
        part: &Part,
    ) {
        let key = mix(sig, code as u64);
        let Some(entry) = Entry::new(code, base, cluster, origin, part) else {
            return;
        };
        let max = self.0.max_bytes;
        if entry.bytes > max / 4 {
            return;
        }
        let mut g = self.lock();
        g.tick += 1;
        let tick = g.tick;
        let bytes = entry.bytes;
        let slots = g.map.entry(key).or_default();
        if !slots.is_empty() {
            // The same cluster computed by two threads at once (or two clusters with one
            // hash, which then share the slot's place in turn): keep what is there.
            return;
        }
        slots.push(Slot {
            entry: Arc::new(entry),
            used: tick,
        });
        g.bytes += bytes;
        if g.bytes > max {
            // Drop the results used least recently, down to three quarters.
            let mut order: Vec<(u64, u64)> = g
                .map
                .iter()
                .flat_map(|(&k, slots)| slots.iter().map(move |s| (s.used, k)))
                .collect();
            order.sort_unstable();
            for (_, k) in order {
                if g.bytes <= max / 4 * 3 {
                    break;
                }
                if let Some(slots) = g.map.remove(&k) {
                    let freed: usize = slots.iter().map(|s| s.entry.bytes).sum();
                    g.bytes = g.bytes.saturating_sub(freed);
                }
            }
        }
    }
}

static VERIFY: AtomicBool = AtomicBool::new(false);

/// Whether reused results are checked (see [`set_memo_verify`]).
pub(crate) fn verify() -> bool {
    VERIFY.load(Ordering::Relaxed)
}

/// Makes every boolean of the process compute each cluster whose result it reuses (from a
/// [`Memo`] or from an equal cluster of the same operation) and panic when the two differ.
/// For tests.
#[doc(hidden)]
pub fn set_memo_verify(yes: bool) {
    VERIFY.store(yes, Ordering::Relaxed);
}
