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
        // Closed fragments as (lo, hi, operand, sign, tag); open ones kept apart.
        let mut f: Vec<(Point, Point, u8, i8, u64)> = Vec::with_capacity(frags.len());
        let mut open_frags = Vec::new();
        for fr in &frags {
            let e = &input[fr.src as usize];
            if e.operand >= 2 {
                open_frags.push(*fr);
            } else if fr.a < fr.b {
                f.push((fr.a, fr.b, e.operand, 1, e.tag));
            } else {
                f.push((fr.b, fr.a, e.operand, -1, e.tag));
            }
        }
        drop(frags);
        f.sort_unstable();
        let mut edges: Vec<MEdge> = Vec::with_capacity(f.len() + open_frags.len());
        let mut i = 0;
        while i < f.len() {
            let (lo, hi) = (f[i].0, f[i].1);
            let mut j = i;
            let mut delta = [0i32; 2];
            while j < f.len() && f[j].0 == lo && f[j].1 == hi {
                delta[f[j].2 as usize] += f[j].3 as i32;
                j += 1;
            }
            if delta != [0, 0] {
                // Tag: first contributor (the group is sorted by operand, sign, tag) whose
                // operand has a non-zero net change.
                let tag = f[i..j]
                    .iter()
                    .find(|x| delta[x.2 as usize] != 0)
                    .map_or(f[i].4, |x| x.4);
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
        for (k, fr) in open_frags.iter().enumerate() {
            let (lo, hi) = if fr.a < fr.b {
                (fr.a, fr.b)
            } else {
                (fr.b, fr.a)
            };
            let tag = input[fr.src as usize].tag;
            edges.push(MEdge {
                lo,
                hi,
                delta: [0, 0],
                tag,
                open: Some(k as u32),
            });
        }
        edges.sort_unstable_by(|a, b| {
            cmp_sweep_edges((a.lo, a.hi), (b.lo, b.hi)).then_with(|| a.open.cmp(&b.open))
        });
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
        Ok(Arrangement {
            edges,
            below,
            open_frags,
        })
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
