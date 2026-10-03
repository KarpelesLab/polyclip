//! Speed comparison of polyclip and Clipper2 (C++, `-O3`) on the spec workloads
//! (POLYGON_LIB.md §5), single thread, coordinates in nanometres.
//!
//! ```text
//! cargo run --release --example timing            # union of 10 000 circles
//! cargo run --release --example timing -- 50000   # the full spec size
//! ```
//!
//! The union workload defaults to 10 000 circles rather than the spec's 50 000 to keep the
//! comparison quick (pass the count as the first argument). Inputs come from the same
//! generator as `examples/perf.rs` (r = 1 mm, 64 vertices, centres in 100 mm x 100 mm).
//!
//! Only the operation itself is timed: for Clipper2, the time inside its `Execute` call
//! (input/output conversion through the C API excluded); for polyclip, the whole public
//! call. Each figure is the median of several runs. Both results are cross-checked by
//! area so a fast wrong answer cannot go unnoticed.

use polyclip::*;
use polyclip_oracle::clipper::{self, Paths};
use polyclip_oracle::compare::{diff, perimeter, rings_of};
use polyclip_oracle::r#gen::{Lcg, circle};
use std::time::{Duration, Instant};

fn median(mut v: Vec<Duration>) -> Duration {
    v.sort();
    v[v.len() / 2]
}

fn runs(n: usize, mut f: impl FnMut() -> Duration) -> Duration {
    median((0..n).map(|_| f()).collect())
}

fn report(name: &str, p: Duration, c: Duration) {
    println!(
        "{name:<44} polyclip {:>10.2?}   Clipper2 {:>10.2?}   ratio {:>5.2}x",
        p,
        c,
        p.as_secs_f64() / c.as_secs_f64()
    );
}

fn check(name: &str, p: &PolygonSet, c: &[Vec<Point>], tol_per_unit: f64) {
    let per = perimeter(rings_of(p).iter().map(|r| r.as_slice()))
        + perimeter(c.iter().map(|r| r.as_slice()));
    let d = diff(p, c, tol_per_unit * per + 4.0);
    println!(
        "{:<44} area polyclip {:.6e}  Clipper2 {:.6e}  xor {:.3e}  {}",
        "",
        d.area_p,
        d.area_c,
        d.xor,
        if d.ok() { "OK" } else { "MISMATCH" }
    );
    assert!(d.ok(), "{name}: results differ beyond tolerance: {d:?}");
}

fn main() {
    let n: usize = std::env::args()
        .nth(1)
        .map(|a| a.parse().expect("circle count"))
        .unwrap_or(10_000);
    let reps = if n > 20_000 { 3 } else { 7 };
    let mut s = Lcg(42);

    // 1. Union of n circles.
    let circles: Vec<Ring> = (0..n)
        .map(|_| {
            let (x, y) = (s.bits() % 100_000_000, s.bits() % 100_000_000);
            circle(x as i64, y as i64, 1_000_000.0, 64)
        })
        .collect();
    let cin = Paths::from_points(circles.iter().map(|r| r.as_slice()));
    let empty = Paths::new();
    let mut pu = PolygonSet::new();
    let tp = runs(reps, || {
        let t = Instant::now();
        pu = union_all(&circles, FillRule::NonZero).unwrap();
        t.elapsed()
    });
    let mut cu = Paths::new();
    let tc = runs(reps, || {
        let (r, dt) = clipper::boolean_timed(Op::Union, &cin, &empty, FillRule::NonZero);
        cu = r;
        dt
    });
    let name = format!("union of {n} circles (64 v)");
    report(&name, tp, tc);
    check(&name, &pu, &cu.to_points(), 2.0);

    // 2. Zone minus 5000 inflated obstacles.
    let side = 100_000_000;
    let zone = Ring::from([(0, 0), (side, 0), (side, side), (0, side)]);
    let obst: Vec<Ring> = (0..5000)
        .map(|_| {
            let (x, y) = (s.bits() % side as u64, s.bits() % side as u64);
            circle(x as i64, y as i64, 300_000.0, 32)
        })
        .collect();
    let zin = Paths::from_points([zone.as_slice()]);
    let oin = Paths::from_points(obst.iter().map(|r| r.as_slice()));
    let mut pd = PolygonSet::new();
    let tp = runs(15, || {
        let t = Instant::now();
        pd = boolean(Op::Difference, &zone, &obst, FillRule::NonZero).unwrap();
        t.elapsed()
    });
    let mut cd = Paths::new();
    let tc = runs(15, || {
        let (r, dt) = clipper::boolean_timed(Op::Difference, &zin, &oin, FillRule::NonZero);
        cd = r;
        dt
    });
    report("zone 100mm - 5000 obstacles (32 v)", tp, tc);
    check("zone", &pd, &cd.to_points(), 2.0);

    // 3. Offset of a 10 000-vertex polygon, round joins, 1 µm arc tolerance.
    let big = circle(0, 0, 50_000_000.0, 10_000);
    let bin = Paths::from_points([big.as_slice()]);
    let delta = 100_000;
    let arc_tol = 1_000;
    let mut po = PolygonSet::new();
    let tp = runs(15, || {
        let t = Instant::now();
        po = offset(&big, delta, Join::Round, ArcTol::new(arc_tol, Side::Inside)).unwrap();
        t.elapsed()
    });
    let mut co = Paths::new();
    let tc = runs(15, || {
        let (r, dt) = clipper::offset_timed(&bin, delta as f64, Join::Round, arc_tol as f64);
        co = r;
        dt
    });
    report("offset 10 000-vertex polygon (+0.1mm round)", tp, tc);
    check("offset", &po, &co.to_points(), arc_tol as f64 + 2.0);
}
