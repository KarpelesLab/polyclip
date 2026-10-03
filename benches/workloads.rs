//! Criterion benchmarks for the spec's performance workloads (POLYGON_LIB.md §5).
//!
//! Coordinates are in nanometres, as in cadlab. Inputs are generated deterministically
//! (fixed-seed LCG) outside the timed closures.
//!
//! * `union_circles/N`: union of N 64-vertex circles of radius 1 mm, heavily overlapping.
//!   The field is scaled with N so every size has the density of the spec case (50 000
//!   circles in 100 mm x 100 mm); target for N = 50 000: < 300 ms.
//! * `zone_minus_obstacles`: one 100 mm x 100 mm zone minus 5 000 inflated obstacles
//!   (32-vertex circles, r = 0.3 mm); target < 50 ms.
//! * `offset_10k`: offset of a 10 000-vertex polygon by +0.1 mm with round joins
//!   (1 µm arc tolerance), convex (circle) and non-convex (wavy star); target < 10 ms.
//!
//! Distance-query benchmarks (`distance_less_than` between two 64-vertex polygons) will
//! be added with the distance API.
//!
//! `cargo bench` runs everything with small sample sizes (a few seconds per benchmark);
//! filter with e.g. `cargo bench --bench workloads -- union_circles/1000`.

use criterion::{
    BenchmarkId, Criterion, SamplingMode, Throughput, criterion_group, criterion_main,
};
use polyclip::*;
use std::f64::consts::PI;
use std::hint::black_box;
use std::time::Duration;

const MM: i64 = 1_000_000;

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    fn below(&mut self, n: i64) -> i64 {
        (self.next() % n as u64) as i64
    }
}

fn circle(cx: i64, cy: i64, r: f64, n: usize) -> Ring {
    (0..n)
        .map(|k| {
            let a = 2.0 * PI * k as f64 / n as f64;
            Point::new(
                cx + (r * a.cos()).round() as i64,
                cy + (r * a.sin()).round() as i64,
            )
        })
        .collect()
}

/// `n` circles with the density of 50 000 circles of r = 1 mm in 100 mm x 100 mm.
fn overlapping_circles(n: usize, seed: u64) -> Vec<Ring> {
    let side = (100.0 * MM as f64 * (n as f64 / 50_000.0).sqrt()) as i64;
    let mut s = Lcg(seed);
    (0..n)
        .map(|_| circle(s.below(side), s.below(side), MM as f64, 64))
        .collect()
}

fn zone_and_obstacles(seed: u64) -> (Ring, Vec<Ring>) {
    let side = 100 * MM;
    let zone = Ring::from([(0, 0), (side, 0), (side, side), (0, side)]);
    let mut s = Lcg(seed);
    let obst = (0..5000)
        .map(|_| circle(s.below(side), s.below(side), 0.3 * MM as f64, 32))
        .collect();
    (zone, obst)
}

/// 10 000-vertex star: radius alternating between 45 mm and 50 mm every 5 vertices.
fn wavy_star(n: usize) -> Ring {
    (0..n)
        .map(|k| {
            let a = 2.0 * PI * k as f64 / n as f64;
            let r = if (k / 5) % 2 == 0 { 50.0 } else { 45.0 } * MM as f64;
            Point::new((r * a.cos()).round() as i64, (r * a.sin()).round() as i64)
        })
        .collect()
}

fn union_circles(c: &mut Criterion) {
    let mut g = c.benchmark_group("union_circles");
    // Long iterations: flat sampling, 10 samples, budget scaled with the size (the 50 000
    // case takes a couple of seconds per iteration today).
    g.sample_size(10);
    g.sampling_mode(SamplingMode::Flat);
    g.warm_up_time(Duration::from_millis(500));
    for n in [1_000usize, 10_000, 50_000] {
        g.measurement_time(Duration::from_secs(match n {
            1_000 => 2,
            10_000 => 6,
            _ => 30,
        }));
        let rings = overlapping_circles(n, 42);
        g.throughput(Throughput::Elements(n as u64));
        g.bench_with_input(BenchmarkId::from_parameter(n), &rings, |b, rings| {
            b.iter(|| union_all(black_box(rings), FillRule::NonZero).unwrap())
        });
    }
    g.finish();
}

fn zone_minus_obstacles(c: &mut Criterion) {
    let (zone, obst) = zone_and_obstacles(7);
    let mut g = c.benchmark_group("zone_minus_obstacles");
    g.sample_size(20);
    g.warm_up_time(Duration::from_millis(500));
    g.measurement_time(Duration::from_secs(3));
    g.bench_function("5000_circles_32v", |b| {
        b.iter(|| {
            boolean(
                Op::Difference,
                black_box(&zone),
                black_box(&obst),
                FillRule::NonZero,
            )
            .unwrap()
        })
    });
    g.finish();
}

fn offset_10k(c: &mut Criterion) {
    let tol = ArcTol::new(1_000, Side::Outside);
    let mut g = c.benchmark_group("offset_10k");
    g.sample_size(20);
    g.warm_up_time(Duration::from_millis(500));
    g.measurement_time(Duration::from_secs(5));
    let convex = circle(0, 0, 50.0 * MM as f64, 10_000);
    let star = wavy_star(10_000);
    for (name, ring) in [("convex", &convex), ("wavy_star", &star)] {
        g.bench_with_input(BenchmarkId::new("round", name), ring, |b, ring| {
            b.iter(|| offset(black_box(ring), MM / 10, Join::Round, tol).unwrap())
        });
    }
    g.finish();
}

criterion_group!(benches, union_circles, zone_minus_obstacles, offset_10k);
criterion_main!(benches);
