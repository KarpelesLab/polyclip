//! Open-path clipping cluster by cluster: only the clip rings near the paths are noded.
//!
//! The locality argument of `cluster.rs` applies to open paths as well: snap rounding is
//! local to segment boxes, so grouping paths and clip rings into clusters (units whose
//! edge boxes come within [`MARGIN`] of each other, transitively) gives clusters whose
//! arrangements are those of the whole input, restricted to them; every clip ring outside a
//! cluster winds a constant number of times around it (its *base winding*, found by
//! locating one point of the cluster). Open edges carry no winding, so only clusters
//! containing paths are computed at all: a few tracks crossing a large zone fill node only
//! the rings they come near.

use crate::arrangement::{Arrangement, InEdge, Noding};
use crate::boolean::FillRule;
use crate::geom::{Point, Rect};
use crate::predicates::orient;
use crate::prepared::Index;

/// Edge boxes within this distance put their units in the same cluster (as in
/// `cluster.rs`).
const MARGIN: i64 = 1;
/// Consecutive edges per chunk box.
const CHUNK: usize = 32;
/// Smaller clip inputs are always computed in one piece.
const MIN_CLIP_EDGES: usize = 1024;

/// A fragment of open edge `src`, from `a` to `b`, and whether it is inside the clip region.
pub(crate) type OpenFrag = (Point, Point, u32, bool);

/// Fragments of the open edges of `edges` (those from `first_open` on) with their inside
/// flags, in the order of [`Arrangement::open_frags`] for the whole input, or `None` when
/// the input does not split usefully (the caller computes it in one piece).
///
/// `ring_starts` are the first edges of the clip rings (all before `first_open`),
/// `path_starts` those of the paths.
pub(crate) fn clustered(
    edges: &[InEdge],
    ring_starts: &[u32],
    first_open: usize,
    path_starts: &[usize],
    rule: FillRule,
) -> Option<Vec<OpenFrag>> {
    if first_open < MIN_CLIP_EDGES || first_open == edges.len() {
        return None;
    }
    // Units: clip rings, then paths (edge ranges, empty ones dropped).
    let mut units: Vec<(u32, u32)> = Vec::new();
    for (i, &s) in ring_starts.iter().enumerate() {
        let e = ring_starts.get(i + 1).map_or(first_open, |&e| e as usize);
        if e > s as usize {
            units.push((s, e as u32));
        }
    }
    let first_path_unit = units.len();
    for (i, &s) in path_starts.iter().enumerate() {
        let e = path_starts.get(i + 1).copied().unwrap_or(edges.len());
        if e > s {
            units.push((s as u32, e as u32));
        }
    }
    if first_path_unit == units.len() {
        return None;
    }
    // Chunks of consecutive edges of a unit, and an index over their boxes.
    let ebox = |k: usize| Rect::new(edges[k].a, edges[k].b);
    let mut chunks: Vec<(Rect, u32, u32, u32)> = Vec::new(); // (box, unit, first, end)
    for (u, &(s, e)) in units.iter().enumerate() {
        let mut k = s;
        while k < e {
            let end = (k + CHUNK as u32).min(e);
            let mut b = ebox(k as usize);
            for j in k + 1..end {
                b = b.union(&ebox(j as usize));
            }
            chunks.push((b, u as u32, k, end));
            k = end;
        }
    }
    let index = Index::of_boxes(chunks.len(), |c| chunks[c].0);
    let mut unit_chunks: Vec<u32> = vec![0; units.len() + 1];
    for c in &chunks {
        unit_chunks[c.1 as usize + 1] += 1;
    }
    for u in 0..units.len() {
        unit_chunks[u + 1] += unit_chunks[u];
    }
    // Clusters grown from every path.
    let mut cluster_of = vec![u32::MAX; units.len()];
    let mut clusters: Vec<Vec<u32>> = Vec::new();
    let budget = first_open * 9 / 10;
    let mut size = 0usize;
    let mut cand: Vec<u32> = Vec::new();
    for p in first_path_unit..units.len() {
        if cluster_of[p] != u32::MAX {
            continue;
        }
        let id = clusters.len() as u32;
        let mut members = vec![p as u32];
        cluster_of[p] = id;
        let mut next = 0;
        while next < members.len() {
            let u = members[next] as usize;
            next += 1;
            if u < first_path_unit {
                size += (units[u].1 - units[u].0) as usize;
                if size > budget {
                    return None;
                }
            }
            for c in unit_chunks[u]..unit_chunks[u + 1] {
                let (cb, _, cs, ce) = chunks[c as usize];
                let g = cb.expand(MARGIN);
                cand.clear();
                index.visit(
                    |b| b.intersects(&g),
                    |o| {
                        let v = chunks[o as usize].1;
                        if cluster_of[v as usize] == u32::MAX && chunks[o as usize].0.intersects(&g)
                        {
                            cand.push(o);
                        }
                    },
                );
                for &o in &cand {
                    let (_, v, os, oe) = chunks[o as usize];
                    if cluster_of[v as usize] != u32::MAX {
                        continue;
                    }
                    let near = (cs..ce).any(|i| {
                        let bi = ebox(i as usize).expand(MARGIN);
                        (os..oe).any(|j| bi.intersects(&ebox(j as usize)))
                    });
                    if near {
                        cluster_of[v as usize] = id;
                        members.push(v);
                    }
                }
            }
        }
        members.sort_unstable();
        clusters.push(members);
    }
    // Base windings: the clip rings outside each cluster around its first path vertex.
    let full = chunks.iter().map(|c| c.0).reduce(|a, b| a.union(&b))?;
    let base: Vec<i32> = clusters
        .iter()
        .enumerate()
        .map(|(id, m)| {
            let q = edges[units
                [*m.iter().find(|&&u| u as usize >= first_path_unit).unwrap() as usize]
                .0 as usize]
                .a;
            let ray = Rect::new(q, Point::new(full.max.x.max(q.x), q.y));
            let mut w = 0i32;
            index.visit(
                |b| b.intersects(&ray),
                |o| {
                    let (b, v, s, e) = chunks[o as usize];
                    if v as usize >= first_path_unit
                        || cluster_of[v as usize] == id as u32
                        || !b.intersects(&ray)
                    {
                        return;
                    }
                    for k in s..e {
                        let (a, b) = (edges[k as usize].a, edges[k as usize].b);
                        if a.y <= q.y {
                            if b.y > q.y && orient(a, b, q) > 0 {
                                w += 1;
                            }
                        } else if b.y <= q.y && orient(a, b, q) < 0 {
                            w -= 1;
                        }
                    }
                },
            );
            w
        })
        .collect();
    // Every cluster on its own.
    let ids: Vec<u32> = (0..clusters.len() as u32).collect();
    let parts: Vec<Vec<OpenFrag>> = crate::par::map_items(&ids, |&id| {
        let m = &clusters[id as usize];
        let b = base[id as usize];
        let mut ce: Vec<InEdge> = Vec::new();
        let mut src: Vec<u32> = Vec::new();
        for &u in m {
            let (s, e) = units[u as usize];
            if (u as usize) >= first_path_unit {
                src.extend(s..e);
            } else {
                ce.extend_from_slice(&edges[s as usize..e as usize]);
            }
        }
        let first = ce.len() as u32;
        ce.extend(src.iter().map(|&k| edges[k as usize]));
        let Ok(arr) = Arrangement::build(&ce, Noding::Snap) else {
            unreachable!("snap rounding never fails")
        };
        let mut inside = vec![false; arr.open_frags.len()];
        for (k, e) in arr.edges.iter().enumerate() {
            if let Some(o) = e.open {
                let (wb, wa) = arr.sides(k);
                inside[o as usize] = rule.is_inside(wb[1] + b) || rule.is_inside(wa[1] + b);
            }
        }
        arr.open_frags
            .iter()
            .zip(inside)
            .map(|(f, ins)| (f.a, f.b, src[(f.src - first) as usize], ins))
            .collect()
    });
    let mut all: Vec<OpenFrag> = parts.into_iter().flatten().collect();
    // Grouped by source edge (stable: fragments of an edge stay in order).
    all.sort_by_key(|f| f.2);
    Some(all)
}
