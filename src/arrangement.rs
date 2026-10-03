//! The noded arrangement of an operation's input edges, with winding numbers.

use crate::geom::Point;
use crate::node::{Crossing, Frag, node_exact, snap_round};
use crate::sweep::{cmp_sweep_edges, sweep};

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
    pub fn build_unwound(input: &[InEdge]) -> Arrangement {
        let segs: Vec<(Point, Point)> = input.iter().map(|e| (e.a, e.b)).collect();
        let frags = snap_round(&segs);
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
