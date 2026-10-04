//! Optional internal parallelism (`rayon` feature). Every helper returns results in the
//! same order as the sequential version, so output never depends on the thread count.

use core::ops::Range;

/// Splits `0..n` into ranges, applies `f` to each (in parallel with the `rayon` feature) and
/// returns the results in range order.
#[cfg(feature = "rayon")]
pub(crate) fn map_ranges<T: Send>(n: usize, f: impl Fn(Range<usize>) -> T + Sync + Send) -> Vec<T> {
    use rayon::prelude::*;
    let threads = rayon::current_num_threads().max(1);
    // Enough chunks for load balancing, few enough to keep per-chunk overhead small.
    let chunks = if n < 4096 {
        1
    } else {
        (threads * 8).min(n / 512).max(1)
    };
    let size = n.div_ceil(chunks).max(1);
    (0..chunks)
        .into_par_iter()
        .map(|k| f((k * size).min(n)..((k + 1) * size).min(n)))
        .collect()
}

/// Sequential fallback: one range.
#[cfg(not(feature = "rayon"))]
pub(crate) fn map_ranges<T: Send>(n: usize, f: impl Fn(Range<usize>) -> T + Sync + Send) -> Vec<T> {
    vec![f(0..n)]
}

/// Applies `f` to every item (in parallel with the `rayon` feature), returning the results
/// in item order.
#[cfg(feature = "rayon")]
pub(crate) fn map_items<T: Sync, U: Send>(
    items: &[T],
    f: impl Fn(&T) -> U + Sync + Send,
) -> Vec<U> {
    use rayon::prelude::*;
    items.par_iter().map(f).collect()
}

/// Sequential fallback.
#[cfg(not(feature = "rayon"))]
pub(crate) fn map_items<T: Sync, U: Send>(
    items: &[T],
    f: impl Fn(&T) -> U + Sync + Send,
) -> Vec<U> {
    items.iter().map(f).collect()
}

/// Sorts by `cmp`, which must order primarily by `x(e)` ascending: elements are first
/// distributed into buckets of nearby `x` (a linear counting pass), then each bucket is
/// sorted (in parallel with the `rayon` feature). Much more cache-friendly than one global
/// comparison sort on large inputs.
pub(crate) fn bucket_sort_by_x<T: Copy + Default + Send + Sync>(
    v: Vec<T>,
    x: impl Fn(&T) -> i64 + Sync,
    cmp: impl Fn(&T, &T) -> core::cmp::Ordering + Sync,
) -> Vec<T> {
    let n = v.len();
    if n < 1 << 14 {
        let mut v = v;
        v.sort_unstable_by(cmp);
        return v;
    }
    let (lo, hi) = v
        .iter()
        .fold((i64::MAX, i64::MIN), |(l, h), e| (l.min(x(e)), h.max(x(e))));
    let span = (hi as i128 - lo as i128 + 1) as u128;
    let nb = (n / 8).max(1);
    let bucket = |e: &T| ((x(e) as i128 - lo as i128) as u128 * nb as u128 / span) as usize;
    let mut start = vec![0usize; nb + 1];
    for e in &v {
        start[bucket(e) + 1] += 1;
    }
    for i in 0..nb {
        start[i + 1] += start[i];
    }
    let mut pos = start.clone();
    let mut out = vec![T::default(); n];
    for e in &v {
        let b = bucket(e);
        out[pos[b]] = *e;
        pos[b] += 1;
    }
    drop(v);
    sort_buckets(&mut out, &start, &cmp);
    out
}

#[cfg(feature = "rayon")]
fn sort_buckets<T: Send>(
    out: &mut [T],
    start: &[usize],
    cmp: &(impl Fn(&T, &T) -> core::cmp::Ordering + Sync),
) {
    use rayon::prelude::*;
    // Split the slice at bucket boundaries into chunks of roughly equal size.
    let n = out.len();
    let chunks = (rayon::current_num_threads() * 8).max(1);
    let mut cuts: Vec<usize> = (1..chunks)
        .map(|k| {
            start[start
                .partition_point(|&s| s < k * n / chunks)
                .min(start.len() - 1)]
        })
        .collect();
    cuts.dedup();
    let mut parts: Vec<(&mut [T], usize)> = Vec::new();
    let mut rest = out;
    let mut base = 0usize;
    for c in cuts {
        if c <= base || c >= n {
            continue;
        }
        let (a, b) = rest.split_at_mut(c - base);
        parts.push((a, base));
        rest = b;
        base = c;
    }
    parts.push((rest, base));
    parts.into_par_iter().for_each(|(part, base)| {
        let end = base + part.len();
        let first = start.partition_point(|&s| s <= base) - 1;
        let mut b = first;
        while b + 1 < start.len() && start[b] < end {
            let (s, e) = (start[b].max(base) - base, start[b + 1].min(end) - base);
            if e > s + 1 {
                part[s..e].sort_unstable_by(cmp);
            }
            b += 1;
        }
    });
}

#[cfg(not(feature = "rayon"))]
fn sort_buckets<T>(out: &mut [T], start: &[usize], cmp: &impl Fn(&T, &T) -> core::cmp::Ordering) {
    for w in start.windows(2) {
        if w[1] > w[0] + 1 {
            out[w[0]..w[1]].sort_unstable_by(cmp);
        }
    }
}
