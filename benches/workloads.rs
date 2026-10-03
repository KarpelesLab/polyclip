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
//! * `distance_less_than`: DRC-style threshold queries between 64-vertex polygons (target
//!   < 1 µs on average including bounding-box rejection), and between near pairs only.
//! * `zone_pipeline`: fracture and triangulation of the zone-fill result.
//! * `pathological`: 50 000 stacked axis-parallel slots and 50 000 diagonal slots (long,
//!   dense, parallel edges).
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

fn distance_queries(c: &mut Criterion) {
    let mut s = Lcg(5);
    let shapes: Vec<Ring> = (0..1000)
        .map(|_| circle(s.below(10 * MM), s.below(10 * MM), 300_000.0, 64))
        .collect();
    let d = 200_000;
    let mut all = Vec::new();
    let mut near = Vec::new();
    for i in 0..shapes.len() {
        let bi = shapes[i].bbox().unwrap();
        for (j, sj) in shapes.iter().enumerate().skip(i + 1).take(100) {
            all.push((i, j));
            if bi.expand(d).intersects(&sj.bbox().unwrap()) {
                near.push((i, j));
            }
        }
    }
    let mut g = c.benchmark_group("distance_less_than");
    g.sample_size(20);
    g.warm_up_time(Duration::from_millis(500));
    g.measurement_time(Duration::from_secs(3));
    for (name, pairs) in [("mixed_pairs", &all), ("near_pairs", &near)] {
        g.throughput(Throughput::Elements(pairs.len() as u64));
        g.bench_with_input(BenchmarkId::new("64v", name), pairs, |b, pairs| {
            b.iter(|| {
                pairs
                    .iter()
                    .filter(|&&(i, j)| {
                        distance_less_than(black_box(&shapes[i]), black_box(&shapes[j]), d)
                    })
                    .count()
            })
        });
    }
    g.finish();
}

fn zone_pipeline(c: &mut Criterion) {
    let (zone, obst) = zone_and_obstacles(7);
    let fill = boolean(Op::Difference, &zone, &obst, FillRule::NonZero).unwrap();
    let mut g = c.benchmark_group("zone_pipeline");
    g.sample_size(10);
    g.warm_up_time(Duration::from_millis(500));
    g.measurement_time(Duration::from_secs(3));
    g.bench_function("fracture", |b| {
        b.iter(|| fracture_set(black_box(&fill)).unwrap())
    });
    g.bench_function("triangulate", |b| {
        b.iter(|| triangulate_set(black_box(&fill)).unwrap())
    });
    g.bench_function("opening_0.1mm", |b| {
        b.iter(|| opening(black_box(&fill), MM / 10, ArcTol::new(1_000, Side::Inside)).unwrap())
    });
    g.finish();
}

fn pathological(c: &mut Criterion) {
    let slots: Vec<Ring> = (0..50_000i64)
        .map(|i| {
            Ring::from([
                (0, 20 * i),
                (1_000_000, 20 * i),
                (1_000_000, 20 * i + 10),
                (0, 20 * i + 10),
            ])
        })
        .collect();
    let diag: Vec<Ring> = (0..50_000i64)
        .map(|i| {
            let (x, y) = (20 * i, -20 * i);
            Ring::from([
                (x, y),
                (x + 700_000, y + 700_000),
                (x + 700_010, y + 699_990),
                (x + 10, y - 10),
            ])
        })
        .collect();
    let mut g = c.benchmark_group("pathological");
    g.sample_size(10);
    g.warm_up_time(Duration::from_millis(500));
    g.measurement_time(Duration::from_secs(3));
    g.bench_function("stacked_slots_50k", |b| {
        b.iter(|| union_all(black_box(&slots), FillRule::NonZero).unwrap())
    });
    g.bench_function("diagonal_slots_50k", |b| {
        b.iter(|| union_all(black_box(&diag), FillRule::NonZero).unwrap())
    });
    g.finish();
}

criterion_group!(
    benches,
    union_circles,
    zone_minus_obstacles,
    offset_10k,
    distance_queries,
    zone_pipeline,
    pathological
);
criterion_main!(benches);
