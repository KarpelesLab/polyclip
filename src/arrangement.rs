//! The noded arrangement of an operation's input edges, with winding numbers.

use crate::assemble::DirEdge;
use crate::geom::Point;
use crate::node::{Crossing, Frag, node_exact, snap_round, snap_round_chunks, snap_round_with};
use crate::sweep::{cmp_sweep_edges, sweep, sweep_band};

/// An input edge. Operands 0 (subject) and 1 (clip) are closed rings; operand 2 marks
/// open-path edges, which take part in noding but carry no winding.
#[derive(Clone, Copy, Debug)]
pub(crate) struct InEdge {
    pub a: Point,
    pub b: Point,
    pub tag: u64,
    pub operand: u8,
}

/// An arrangement edge, `lo < hi`.
#[derive(Clone, Copy, Debug)]
pub(crate) struct MEdge {
    pub lo: Point,
    pub hi: Point,
    /// Winding change for each closed operand when crossing from below to above.
    pub delta: [i32; 2],
    pub tag: u64,
    /// For open-path fragments: index into [`Arrangement::open_frags`].
    pub open: Option<u32>,
}

pub(crate) enum Noding {
    /// Snap rounding (always succeeds).
    Snap,
    /// Exact noding (fails on a proper crossing).
    Exact,
}

pub(crate) struct Arrangement {
    /// Edges in sweep order. Closed edges are unique; open fragments follow any coincident
    /// closed edge.
    pub edges: Vec<MEdge>,
    /// Winding numbers just below each edge.
    pub below: Vec<[i32; 2]>,
    /// Open-path fragments, grouped by source edge in input order.
    pub open_frags: Vec<Frag>,
}

impl Arrangement {
    pub fn build(input: &[InEdge], noding: Noding) -> Result<Arrangement, Crossing> {
        let segs: Vec<(Point, Point)> = input.iter().map(|e| (e.a, e.b)).collect();
        let frags = match noding {
            Noding::Snap => snap_round(&segs),
            Noding::Exact => node_exact(&segs)?,
        };
        drop(segs);
        Ok(Self::from_frags(input, frags))
    }

    /// Like [`build`](Self::build) with snap rounding, but without computing the winding
    /// numbers (`below` stays empty): for callers that run their own sweep.
    pub fn build_unwound(input: &[InEdge], engine: Engine) -> Arrangement {
        let segs: Vec<(Point, Point)> = input.iter().map(|e| (e.a, e.b)).collect();
        let frags = snap_round_with(&segs, engine);
        drop(segs);
        Self::merge(input, frags)
    }

    /// Builds the arrangement from already-noded fragments of `input`.
    pub fn from_frags(input: &[InEdge], frags: Vec<Frag>) -> Arrangement {
        let mut arr = Self::merge(input, frags);
        arr.wind();
        arr
    }

    /// Merges coincident fragments into sorted arrangement edges (no windings yet).
    fn merge(input: &[InEdge], frags: Vec<Frag>) -> Arrangement {
        #[derive(Clone, Copy, Default)]
        struct F {
            lo: Point,
            hi: Point,
            tag: u64,
            /// `u32::MAX` for closed fragments, else the open fragment index.
            open: u32,
            operand: u8,
            sign: i8,
        }
        let mut f: Vec<F> = Vec::with_capacity(frags.len());
        let mut open_frags = Vec::new();
        for fr in &frags {
            let e = &input[fr.src as usize];
            let (lo, hi, sign) = if fr.a < fr.b {
                (fr.a, fr.b, 1)
            } else {
                (fr.b, fr.a, -1)
            };
            let open = if e.operand >= 2 {
                open_frags.push(*fr);
                open_frags.len() as u32 - 1
            } else {
                u32::MAX
            };
            f.push(F {
                lo,
                hi,
                tag: e.tag,
                open,
                operand: e.operand,
                sign,
            });
        }
        drop(frags);
        // One sort into sweep order. Coincident fragments are adjacent (same start, same
        // direction); closed ones come first, ordered by operand, sign and tag so the merge
        // below is deterministic.
        let f = crate::par::bucket_sort_by_x(
            f,
            |e| e.lo.x,
            |a, b| {
                cmp_sweep_edges((a.lo, a.hi), (b.lo, b.hi))
                    .then_with(|| a.hi.cmp(&b.hi))
                    .then_with(|| (b.open == u32::MAX).cmp(&(a.open == u32::MAX)))
                    .then_with(|| {
                        (a.operand, a.sign, a.tag, a.open).cmp(&(b.operand, b.sign, b.tag, b.open))
                    })
            },
        );
        let mut edges: Vec<MEdge> = Vec::with_capacity(f.len());
        let mut i = 0;
        while i < f.len() {
            let (lo, hi) = (f[i].lo, f[i].hi);
            if f[i].open != u32::MAX {
                edges.push(MEdge {
                    lo,
                    hi,
                    delta: [0, 0],
                    tag: f[i].tag,
                    open: Some(f[i].open),
                });
                i += 1;
                continue;
            }
            let mut j = i;
            let mut delta = [0i32; 2];
            while j < f.len() && f[j].open == u32::MAX && f[j].lo == lo && f[j].hi == hi {
                delta[f[j].operand as usize] += f[j].sign as i32;
                j += 1;
            }
            if delta != [0, 0] {
                // Tag: first contributor (by operand, sign, tag) whose operand has a non-zero
                // net change.
                let tag = f[i..j]
                    .iter()
                    .find(|x| delta[x.operand as usize] != 0)
                    .map_or(f[i].tag, |x| x.tag);
                edges.push(MEdge {
                    lo,
                    hi,
                    delta,
                    tag,
                    open: None,
                });
            }
            i = j;
        }
        drop(f);
        Arrangement {
            edges,
            below: Vec::new(),
            open_frags,
        }
    }

    /// Computes the winding numbers below every edge.
    fn wind(&mut self) {
        let edges = &self.edges;
        let segs: Vec<(Point, Point)> = edges.iter().map(|e| (e.lo, e.hi)).collect();
        let mut below = vec![[0i32; 2]; edges.len()];
        sweep(&segs, |e, b| {
            if let Some(b) = b {
                let b = b as usize;
                below[e as usize] = [
                    below[b][0] + edges[b].delta[0],
                    below[b][1] + edges[b].delta[1],
                ];
            }
        });
        self.below = below;
    }

    /// Winding numbers below and above edge `k`. For an open fragment lying on a closed
    /// edge, those of the closed edge (so both adjacent faces are seen).
    pub fn sides(&self, k: usize) -> ([i32; 2], [i32; 2]) {
        let e = &self.edges[k];
        let mut c = k;
        if e.open.is_some() {
            while c > 0 {
                let p = &self.edges[c - 1];
                if p.lo != e.lo || p.hi != e.hi {
                    break;
                }
                c -= 1;
                if p.open.is_none() {
                    break;
                }
            }
            if self.edges[c].open.is_some() {
                c = k;
            }
        }
        let wb = self.below[c];
        let d = self.edges[c].delta;
        (wb, [wb[0] + d[0], wb[1] + d[1]])
    }
}

/// How booleans compute their arrangement and boundary. The output is the same in every
/// mode; for tests and benchmarks.
#[doc(hidden)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Engine {
    /// Parallel phases sized for the machine (the default).
    #[default]
    Auto,
    /// The original sequential implementation, kept for comparison.
    Reference,
    /// The parallel phases split into as many pieces as possible, so that the code merging
    /// them is exercised on small inputs too.
    Split,
}

std::thread_local! {
    static ENGINE: core::cell::Cell<Engine> = const { core::cell::Cell::new(Engine::Auto) };
}

/// Selects how booleans started on the current thread compute their result (they then run
/// in one piece, without splitting into clusters). The output is the same either way; for
/// tests and benchmarks.
#[doc(hidden)]
pub fn set_engine(e: Engine) {
    ENGINE.with(|m| m.set(e));
}

/// The engine selected on the current thread.
pub(crate) fn engine() -> Engine {
    ENGINE.with(|m| m.get())
}

/// Minimal union-find over region ids.
struct Regions {
    parent: Vec<u32>,
}

impl Regions {
    fn add(&mut self) -> u32 {
        self.parent.push(self.parent.len() as u32);
        self.parent.len() as u32 - 1
    }
    fn find(&mut self, mut x: u32) -> u32 {
        while self.parent[x as usize] != x {
            let p = self.parent[self.parent[x as usize] as usize];
            self.parent[x as usize] = p;
            x = p;
        }
        x
    }
    fn union(&mut self, a: u32, b: u32) {
        let (a, b) = (self.find(a), self.find(b));
        if a != b {
            // Keep the smaller id as root (deterministic).
            let (lo, hi) = if a < b { (a, b) } else { (b, a) };
            self.parent[hi as usize] = lo;
        }
    }
}

/// Computes the directed boundary edges of the region `{ inside(winding) }`, in sweep
/// order, each with the connected regions on its two sides.
///
/// One sweep computes winding numbers and, with a union-find over the gaps between
/// consecutive status edges, the connected components of the result and of its complement:
/// gaps on both sides of a non-boundary edge belong to the same component, and gaps meeting
/// at a vertex where edges only end merge. Region 0 is the unbounded exterior; other region
/// ids only identify components (their values carry no meaning).
///
/// With several threads the sweep is split into vertical bands swept independently, each
/// starting from the status at its left side (whose windings are prefix sums); the bands'
/// regions are then joined through the gaps they share at their common sides.
pub(crate) fn boundary(
    input: &[InEdge],
    inside: impl Fn([i32; 2]) -> bool + Sync,
    engine: Engine,
) -> Vec<DirEdge> {
    if engine == Engine::Reference || input.iter().any(|e| e.operand >= 2) {
        return boundary_reference(input, inside);
    }
    let segs: Vec<(Point, Point)> = crate::par::map_slice(input, |e| (e.a, e.b));
    let frags = snap_round_chunks(&segs, engine);
    drop(segs);
    let lean = merge_lean(input, frags, engine);
    let mut his: Vec<Point> = crate::par::map_slice(&lean.segs, |e| e.1);
    crate::par::sort_unstable(&mut his);
    let cuts = band_cuts(&lean.segs, &his, engine);
    if cuts.len() <= 2 {
        return boundary_sequential(&lean, &his, &inside);
    }
    boundary_banded(&lean, &his, &cuts, &inside)
}

/// The original one-sweep implementation of [`boundary`].
fn boundary_reference(input: &[InEdge], inside: impl Fn([i32; 2]) -> bool) -> Vec<DirEdge> {
    let arr = Arrangement::build_unwound(input, Engine::Reference);
    let lean = Lean {
        segs: arr.edges.iter().map(|e| (e.lo, e.hi)).collect(),
        delta: arr.edges.iter().map(|e| e.delta).collect(),
        tag: arr.edges.iter().map(|e| e.tag).collect(),
    };
    let mut his: Vec<Point> = lean.segs.iter().map(|e| e.1).collect();
    his.sort_unstable();
    boundary_sequential(&lean, &his, &inside)
}

/// The arrangement edges of closed rings in sweep order, as [`Arrangement::build_unwound`]
/// computes them (open paths aside), in the compact form [`boundary`] needs.
struct Lean {
    segs: Vec<(Point, Point)>,
    delta: Vec<[i32; 2]>,
    tag: Vec<u64>,
}

/// [`Arrangement::merge`] for closed edges only, from fragments in chunks: the fragments
/// are distributed into buckets of `lo.x` ranges, and every bucket is sorted and merged on
/// its own (in parallel with the `rayon` feature). Coincident fragments share their `lo`, so
/// they land in the same bucket and the result is the same as with one global sort.
fn merge_lean(input: &[InEdge], chunks: Vec<Vec<Frag>>, engine: Engine) -> Lean {
    /// A fragment oriented `lo -> hi`, with its source edge.
    #[derive(Clone, Copy, Default)]
    struct F {
        lo: Point,
        hi: Point,
        src: u32,
        sign: i8,
    }
    let total: usize = chunks.iter().map(|c| c.len()).sum();
    let want = match engine {
        Engine::Split => 7,
        // Several buckets even on one thread: smaller sorts, less memory at once.
        _ if total >= 1 << 16 => (crate::par::threads() * 8).max(16),
        _ => 1,
    };
    // Bucket bounds: quantiles of a sample of the `lo.x` values.
    let lo_x = |f: &Frag| f.a.x.min(f.b.x);
    let mut bounds: Vec<i64> = Vec::new();
    if want > 1 {
        let step = (total / (want * 64)).max(1);
        let mut xs: Vec<i64> = chunks.iter().flatten().step_by(step).map(lo_x).collect();
        xs.sort_unstable();
        for k in 1..want {
            if let Some(&x) = xs.get(k * xs.len() / want) {
                bounds.push(x);
            }
        }
        bounds.dedup();
    }
    let nbk = bounds.len() + 1;
    // Distribute (each chunk of fragments is freed once converted).
    let parts: Vec<Vec<Vec<F>>> = crate::par::map_vec(chunks, |ch| {
        let mut v: Vec<Vec<F>> = vec![Vec::new(); nbk];
        if nbk == 1 {
            v[0].reserve_exact(ch.len());
        }
        for fr in ch {
            let (lo, hi, sign) = if fr.a < fr.b {
                (fr.a, fr.b, 1)
            } else {
                (fr.b, fr.a, -1)
            };
            let b = bounds.partition_point(|&x| x <= lo.x);
            v[b].push(F {
                lo,
                hi,
                src: fr.src,
                sign,
            });
        }
        v
    });
    let mut buckets: Vec<Vec<Vec<F>>> = (0..nbk).map(|_| Vec::new()).collect();
    for p in parts {
        for (b, v) in p.into_iter().enumerate() {
            buckets[b].push(v);
        }
    }
    let key = |f: &F| {
        let e = &input[f.src as usize];
        (e.operand, f.sign, e.tag)
    };
    let outs: Vec<Lean> = crate::par::map_vec(buckets, |lists| {
        let mut f: Vec<F> = if lists.len() == 1 {
            lists.into_iter().next().unwrap_or_default()
        } else {
            let mut f = Vec::with_capacity(lists.iter().map(|l| l.len()).sum());
            for l in lists {
                f.extend_from_slice(&l);
            }
            f
        };
        // Sweep order; coincident fragments are adjacent, ordered by operand, sign and tag
        // so the merge below is deterministic.
        f = crate::par::bucket_sort_by_x(
            f,
            |e| e.lo.x,
            |a, b| {
                cmp_sweep_edges((a.lo, a.hi), (b.lo, b.hi))
                    .then_with(|| a.hi.cmp(&b.hi))
                    .then_with(|| key(a).cmp(&key(b)))
            },
        );
        let mut out = Lean {
            segs: Vec::with_capacity(f.len()),
            delta: Vec::with_capacity(f.len()),
            tag: Vec::with_capacity(f.len()),
        };
        let mut i = 0;
        while i < f.len() {
            let (lo, hi) = (f[i].lo, f[i].hi);
            let mut j = i;
            let mut delta = [0i32; 2];
            while j < f.len() && f[j].lo == lo && f[j].hi == hi {
                delta[input[f[j].src as usize].operand as usize] += f[j].sign as i32;
                j += 1;
            }
            if delta != [0, 0] {
                // Tag: first contributor (by operand, sign, tag) whose operand has a non-zero
                // net change.
                let first = f[i..j]
                    .iter()
                    .find(|x| delta[input[x.src as usize].operand as usize] != 0)
                    .unwrap_or(&f[i]);
                out.segs.push((lo, hi));
                out.delta.push(delta);
                out.tag.push(input[first.src as usize].tag);
            }
            i = j;
        }
        out
    });
    let mut segs = Vec::with_capacity(outs.len());
    let mut delta = Vec::with_capacity(outs.len());
    let mut tag = Vec::with_capacity(outs.len());
    for o in outs {
        segs.push(o.segs);
        delta.push(o.delta);
        tag.push(o.tag);
    }
    Lean {
        segs: crate::par::concat_vecs(segs),
        delta: crate::par::concat_vecs(delta),
        tag: crate::par::concat_vecs(tag),
    }
}

fn boundary_sequential(
    lean: &Lean,
    his: &[Point],
    inside: &impl Fn([i32; 2]) -> bool,
) -> Vec<DirEdge> {
    let segs = &lean.segs;
    let n = segs.len();
    let mut below_w = vec![[0i32; 2]; n];
    let mut gap_above = vec![0u32; n];
    let mut reg = Regions { parent: vec![0] };
    let mut out: Vec<DirEdge> = Vec::new();
    let delta = |k: usize| lean.delta[k];
    sweep_band(segs, his, 0..n, &[], |below, _, ending, starting| {
        let g_below = below.map_or(0, |b| gap_above[b as usize]);
        let g_above = ending.last().map_or(g_below, |&e| gap_above[e as usize]);
        if starting.is_empty() {
            if !ending.is_empty() {
                reg.union(g_below, g_above);
            }
            return;
        }
        let mut gb = g_below;
        let mut wb = below.map_or([0, 0], |b| {
            let b = b as usize;
            [below_w[b][0] + delta(b)[0], below_w[b][1] + delta(b)[1]]
        });
        let last = starting.end - 1;
        for k in starting {
            let ku = k as usize;
            let d = delta(ku);
            let wa = [wb[0] + d[0], wb[1] + d[1]];
            let ga = if k == last { g_above } else { reg.add() };
            below_w[ku] = wb;
            gap_above[ku] = ga;
            let (ib, ia) = (inside(wb), inside(wa));
            if ib == ia {
                reg.union(gb, ga);
            } else {
                let (lo, hi) = segs[ku];
                // Interior on the left: above for lo->hi.
                let (from, to) = if ia { (lo, hi) } else { (hi, lo) };
                out.push(DirEdge {
                    from,
                    to,
                    tag: lean.tag[ku],
                    below: gb,
                    above: ga,
                });
            }
            gb = ga;
            wb = wa;
        }
    });
    for e in out.iter_mut() {
        e.below = reg.find(e.below);
        e.above = reg.find(e.above);
    }
    out
}

/// Band boundaries for [`boundary_banded`] as edge indices `0 = c_0 < c_1 < ... = n`; band
/// `k` holds the vertices `v` with `lo(c_k) <= v < lo(c_{k+1})`. A single band means one
/// sweep.
fn band_cuts(segs: &[(Point, Point)], his: &[Point], engine: Engine) -> Vec<usize> {
    let n = segs.len();
    let want = match engine {
        Engine::Reference => 1,
        Engine::Split => n / 2,
        Engine::Auto if n >= 1 << 15 && crate::par::threads() > 1 => crate::par::threads() * 4,
        Engine::Auto => 1,
    };
    let mut cuts = vec![0usize];
    if want > 1 {
        for k in 1..want {
            // Start of a vertex's group of edges.
            let mut c = k * n / want;
            while c > 0 && c < n && segs[c].0 == segs[c - 1].0 {
                c += 1;
            }
            if c >= n || c <= *cuts.last().unwrap_or(&0) {
                continue;
            }
            // Edges in the status at the cut (lo < v <= hi), which the band starting there
            // must sort: keep that small relative to the band.
            let width = c - his.partition_point(|h| *h < segs[c].0);
            if engine == Engine::Auto && width * want * 4 > n {
                continue;
            }
            cuts.push(c);
        }
    }
    cuts.push(n);
    cuts
}

/// What a band of [`boundary_banded`] reports.
struct BandOut {
    /// Boundary edges, regions numbered locally (0 is the exterior, `1..=init.len()` the
    /// gaps above the initial status edges).
    out: Vec<DirEdge>,
    /// Local union-find (no entry above its index).
    parent: Vec<u32>,
    /// The local gap above every edge of the final status (bottom to top).
    fin_gap: Vec<u32>,
}

fn boundary_banded(
    lean: &Lean,
    his: &[Point],
    cuts: &[usize],
    inside: &(impl Fn([i32; 2]) -> bool + Sync),
) -> Vec<DirEdge> {
    let segs = &lean.segs;
    let nb = cuts.len() - 1;
    let n = segs.len();
    let vtx = |k: usize| segs[cuts[k]].0;
    // Initial status of every band: the edges with lo < v_k <= hi.
    let pairs: Vec<Vec<(u32, u32)>> = crate::par::map_ranges(n, |range| {
        let mut out = Vec::new();
        let mut b = cuts.partition_point(|&c| c <= range.start).max(1) - 1;
        for e in range {
            while cuts[b + 1] <= e {
                b += 1;
            }
            let hi = segs[e].1;
            let mut k = b + 1;
            while k < nb && vtx(k) <= hi {
                out.push((k as u32, e as u32));
                k += 1;
            }
        }
        out
    });
    let mut istart = vec![0usize; nb + 1];
    for &(k, _) in pairs.iter().flatten() {
        istart[k as usize + 1] += 1;
    }
    for k in 0..nb {
        istart[k + 1] += istart[k];
    }
    let mut init: Vec<u32> = vec![0; istart[nb]];
    let mut pos = istart.clone();
    for &(k, e) in pairs.iter().flatten() {
        init[pos[k as usize]] = e;
        pos[k as usize] += 1;
    }
    drop(pairs);
    let ids: Vec<usize> = (0..nb).collect();
    let bands: Vec<BandOut> = crate::par::map_items(&ids, |&k| {
        let mut status: Vec<u32> = init[istart[k]..istart[k + 1]].to_vec();
        status.sort_unstable_by(|&a, &b| crate::sweep::cmp_status(segs, a, b));
        let h0 = if k == 0 {
            0
        } else {
            his.partition_point(|h| *h < vtx(k))
        };
        let h1 = if k + 1 == nb {
            his.len()
        } else {
            his.partition_point(|h| *h < vtx(k + 1))
        };
        sweep_one_band(lean, &his[h0..h1], cuts[k]..cuts[k + 1], &status, inside)
    });
    // Global regions: band k's local id l > 0 becomes off[k] + l.
    let mut off = Vec::with_capacity(nb);
    let mut total = 1usize;
    for b in &bands {
        off.push(total as u32 - 1);
        total += b.parent.len() - 1;
    }
    let g = |k: usize, l: u32| if l == 0 { 0 } else { off[k] + l };
    let mut reg = Regions {
        parent: vec![0; total],
    };
    for (k, b) in bands.iter().enumerate() {
        for (l, &p) in b.parent.iter().enumerate().skip(1) {
            reg.parent[g(k, l as u32) as usize] = g(k, p);
        }
    }
    // Join the gaps the bands share: band k - 1's final status is band k's initial one.
    for k in 1..nb {
        let prev = &bands[k - 1];
        debug_assert_eq!(prev.fin_gap.len(), istart[k + 1] - istart[k]);
        for (i, &gap) in prev.fin_gap.iter().enumerate() {
            reg.union(g(k - 1, gap), g(k, i as u32 + 1));
        }
    }
    let mut out = Vec::with_capacity(bands.iter().map(|b| b.out.len()).sum());
    for (k, b) in bands.into_iter().enumerate() {
        for mut e in b.out {
            e.below = reg.find(g(k, e.below));
            e.above = reg.find(g(k, e.above));
            out.push(e);
        }
    }
    out
}

/// Sweeps one band of [`boundary_banded`]: `starting` are the edges starting in it, `his`
/// the sorted endpoints ending in it and `init` the status at its left side.
fn sweep_one_band(
    lean: &Lean,
    his: &[Point],
    starting: core::ops::Range<usize>,
    init: &[u32],
    inside: &impl Fn([i32; 2]) -> bool,
) -> BandOut {
    let segs = &lean.segs;
    let c0 = starting.start;
    let delta = |k: usize| lean.delta[k];
    // Windings below and gaps above: the band's own edges by index, initial ones by slot.
    let mut below_w = vec![[0i32; 2]; starting.len()];
    let mut gap_above = vec![0u32; starting.len()];
    let mut init_w = Vec::with_capacity(init.len());
    let mut w = [0i32; 2];
    for &e in init {
        init_w.push(w);
        let d = delta(e as usize);
        w = [w[0] + d[0], w[1] + d[1]];
    }
    let mut slot: Vec<(u32, u32)> = init
        .iter()
        .enumerate()
        .map(|(i, &e)| (e, i as u32))
        .collect();
    slot.sort_unstable();
    let find_slot = |e: u32| slot[slot.partition_point(|x| x.0 < e)].1 as usize;
    // Gap above an edge: the initial ones are 1..=init.len().
    let gap = |e: u32, gap_above: &[u32]| {
        if e as usize >= c0 {
            gap_above[e as usize - c0]
        } else {
            find_slot(e) as u32 + 1
        }
    };
    let mut reg = Regions {
        parent: (0..=init.len() as u32).collect(),
    };
    let mut out: Vec<DirEdge> = Vec::new();
    let fin = sweep_band(segs, his, starting, init, |below, _, ending, starting| {
        let g_below = below.map_or(0, |b| gap(b, &gap_above));
        let g_above = ending.last().map_or(g_below, |&e| gap(e, &gap_above));
        if starting.is_empty() {
            if !ending.is_empty() {
                reg.union(g_below, g_above);
            }
            return;
        }
        let mut gb = g_below;
        let mut wb = below.map_or([0, 0], |b| {
            let wbb = if b as usize >= c0 {
                below_w[b as usize - c0]
            } else {
                init_w[find_slot(b)]
            };
            let d = delta(b as usize);
            [wbb[0] + d[0], wbb[1] + d[1]]
        });
        let last = starting.end - 1;
        for k in starting {
            let ku = k as usize;
            let d = delta(ku);
            let wa = [wb[0] + d[0], wb[1] + d[1]];
            let ga = if k == last { g_above } else { reg.add() };
            below_w[ku - c0] = wb;
            gap_above[ku - c0] = ga;
            let (ib, ia) = (inside(wb), inside(wa));
            if ib == ia {
                reg.union(gb, ga);
            } else {
                let (lo, hi) = segs[ku];
                let (from, to) = if ia { (lo, hi) } else { (hi, lo) };
                out.push(DirEdge {
                    from,
                    to,
                    tag: lean.tag[ku],
                    below: gb,
                    above: ga,
                });
            }
            gb = ga;
            wb = wa;
        }
    });
    let fin_gap = fin.iter().map(|&e| gap(e, &gap_above)).collect();
    BandOut {
        out,
        parent: reg.parent,
        fin_gap,
    }
}
