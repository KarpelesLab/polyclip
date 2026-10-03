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

/// Unstable sort (parallel with the `rayon` feature). Callers use comparators under which
/// equal elements are indistinguishable, so the result is the same either way.
#[cfg(feature = "rayon")]
pub(crate) fn sort_unstable_by<T: Send>(
    v: &mut [T],
    cmp: impl Fn(&T, &T) -> core::cmp::Ordering + Sync,
) {
    use rayon::prelude::*;
    if v.len() > 1 << 15 {
        v.par_sort_unstable_by(cmp);
    } else {
        v.sort_unstable_by(cmp);
    }
}

/// Sequential fallback.
#[cfg(not(feature = "rayon"))]
pub(crate) fn sort_unstable_by<T: Send>(
    v: &mut [T],
    cmp: impl Fn(&T, &T) -> core::cmp::Ordering + Sync,
) {
    v.sort_unstable_by(cmp);
}
