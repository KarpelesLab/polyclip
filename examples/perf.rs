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
    if std::env::args().any(|a| a == "wide") {
        if std::env::var("PITCH").is_ok() {
            return diag_only();
        }
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
        let t = Instant::now();
        let u = union_all(&slots, FillRule::NonZero).unwrap();
        println!(
            "union of 50k stacked slots: {:?} -> {}",
            t.elapsed(),
            u.len()
        );
        // The same, rotated by 45 degrees (no axis separates the segments).
        let pitch: i64 = std::env::var("PITCH")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(20);
        let diag: Vec<Ring> = (0..50_000i64)
            .map(|i| {
                let (x, y) = (pitch * i, -pitch * i);
                Ring::from([
                    (x, y),
                    (x + 700_000, y + 700_000),
                    (x + 700_010, y + 699_990),
                    (x + 10, y - 10),
                ])
            })
            .collect();
        let t = Instant::now();
        let u = union_all(&diag, FillRule::NonZero).unwrap();
        println!(
            "union of 50k diagonal slots: {:?} -> {}",
            t.elapsed(),
            u.len()
        );
        return;
    }
    if std::env::args().any(|a| a == "zone") {
        let mut s = 42u64;
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
        for _ in 0..20 {
            std::hint::black_box(boolean(Op::Difference, &zone, &obst, FillRule::NonZero).unwrap());
        }
        println!("zone x20: {:?}", t.elapsed() / 20);
        return;
    }
    if std::env::args().any(|a| a == "dist") {
        distance_bench();
        return;
    }
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
    let t = Instant::now();
    let f = fracture_set(&d).unwrap();
    println!(
        "fracture zone: {:?} -> {} verts",
        t.elapsed(),
        f.iter().map(|r| r.len()).sum::<usize>()
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

fn diag_only() {
    let pitch: i64 = std::env::var("PITCH").unwrap().parse().unwrap();
    let diag: Vec<Ring> = (0..20_000i64)
        .map(|i| {
            let (x, y) = (pitch * i, -pitch * i);
            Ring::from([
                (x, y),
                (x + 700_000, y + 700_000),
                (x + 700_100, y + 699_900),
                (x + 100, y - 100),
            ])
        })
        .collect();
    let t = Instant::now();
    let u = union_all(&diag, FillRule::NonZero).unwrap();
    println!(
        "union of 20k diagonal tracks, pitch {pitch}: {:?} -> {}",
        t.elapsed(),
        u.len()
    );
}

#[allow(dead_code)]
fn distance_bench() {
    let mut s = 7u64;
    let shapes: Vec<Ring> = (0..2000)
        .map(|_| {
            circle(
                (lcg(&mut s) % 10_000_000) as i64,
                (lcg(&mut s) % 10_000_000) as i64,
                300_000.0,
                64,
            )
        })
        .collect();
    let t = Instant::now();
    let mut hits = 0;
    let mut n = 0;
    for i in 0..shapes.len() {
        for j in (i + 1)..shapes.len().min(i + 200) {
            n += 1;
            if distance_less_than(&shapes[i], &shapes[j], 200_000) {
                hits += 1;
            }
        }
    }
    let e = t.elapsed();
    println!(
        "distance_less_than: {} queries, {} hits, {:?}/query",
        n,
        hits,
        e / n as u32
    );
    // Near pairs only (bboxes overlap).
    let boxes: Vec<Rect> = shapes.iter().map(|s| s.bbox().unwrap()).collect();
    let mut pairs = Vec::new();
    for i in 0..shapes.len() {
        let a = boxes[i].expand(200_000);
        for (j, b) in boxes.iter().enumerate() {
            if i != j && a.intersects(b) {
                pairs.push((i, j));
            }
        }
    }
    let t = Instant::now();
    for &(i, j) in &pairs {
        std::hint::black_box(distance_less_than(&shapes[i], &shapes[j], 200_000));
    }
    println!(
        "distance_less_than near pairs: {} queries, {:?}/query",
        pairs.len(),
        t.elapsed() / pairs.len().max(1) as u32
    );
}
