use polyclip::*;
use std::time::Instant;

fn lcg(s: &mut u64) -> u64 {
    *s = s
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    *s >> 33
}

fn circle(cx: i64, cy: i64, r: f64, n: usize) -> Ring {
    (0..n)
        .map(|k| {
            let a = 2.0 * std::f64::consts::PI * k as f64 / n as f64;
            Point::new(
                cx + (r * a.cos()).round() as i64,
                cy + (r * a.sin()).round() as i64,
            )
        })
        .collect()
}

fn main() {
    let mut s = 42u64;
    let n: usize = std::env::args()
        .nth(1)
        .map(|a| a.parse().unwrap())
        .unwrap_or(50_000);
    let circles: Vec<Ring> = (0..n)
        .map(|_| {
            circle(
                (lcg(&mut s) % 100_000_000) as i64,
                (lcg(&mut s) % 100_000_000) as i64,
                1_000_000.0,
                64,
            )
        })
        .collect();
    let t = Instant::now();
    let u = union_all(&circles, FillRule::NonZero).unwrap();
    println!(
        "union {} circles: {:?} -> {} polys, {} verts",
        n,
        t.elapsed(),
        u.len(),
        u.iter().map(|p| p.vertex_count()).sum::<usize>()
    );

    let zone = Ring::from([
        (0, 0),
        (100_000_000, 0),
        (100_000_000, 100_000_000),
        (0, 100_000_000),
    ]);
    let obst: Vec<Ring> = (0..5000)
        .map(|_| {
            circle(
                (lcg(&mut s) % 100_000_000) as i64,
                (lcg(&mut s) % 100_000_000) as i64,
                300_000.0,
                32,
            )
        })
        .collect();
    let t = Instant::now();
    let d = boolean(Op::Difference, &zone, &obst, FillRule::NonZero).unwrap();
    println!(
        "zone - 5000 obstacles: {:?} -> {} polys",
        t.elapsed(),
        d.len()
    );

    let big = circle(0, 0, 50_000_000.0, 10_000);
    let t = Instant::now();
    let o = offset(&big, 100_000, Join::Round, ArcTol::new(1000, Side::Outside)).unwrap();
    println!(
        "offset 10k-vertex polygon: {:?} -> {} verts",
        t.elapsed(),
        o[0].outer.len()
    );
}
