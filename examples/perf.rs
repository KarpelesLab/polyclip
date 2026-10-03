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

fn checksum(ps: &PolygonSet) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut h = std::collections::hash_map::DefaultHasher::new();
    ps.hash(&mut h);
    h.finish()
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
    if std::env::args().any(|a| a == "offset") {
        let big = circle(0, 0, 50_000_000.0, 10_000);
        // A non-convex star as well.
        let star: Ring = (0..10_000)
            .map(|k| {
                let a = 2.0 * std::f64::consts::PI * k as f64 / 10_000.0;
                let r = if k % 2 == 0 {
                    50_000_000.0
                } else {
                    45_000_000.0
                };
                Point::new((r * a.cos()).round() as i64, (r * a.sin()).round() as i64)
            })
            .collect();
        let only_circle = std::env::args().any(|a| a == "circle");
        for (name, p) in [("circle", &big), ("star", &star)] {
            if only_circle && name == "star" {
                continue;
            }
            let t = Instant::now();
            for _ in 0..20 {
                std::hint::black_box(
                    offset(p, 100_000, Join::Round, ArcTol::new(1000, Side::Outside)).unwrap(),
                );
            }
            println!("offset {name} x20: {:?}", t.elapsed() / 20);
            let t = Instant::now();
            for _ in 0..20 {
                std::hint::black_box(union_all(p, FillRule::NonZero).unwrap());
            }
            println!("  normalize only: {:?}", t.elapsed() / 20);
        }
        return;
    }
    if std::env::args().any(|a| a == "circles") {
        let mut s = 42u64;
        let circles: Vec<Ring> = (0..50_000)
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
            "union 50k circles: {:?} -> {} (checksum {:x})",
            t.elapsed(),
            u.len(),
            checksum(&u)
        );
        return;
    }
    if std::env::args().any(|a| a == "opening") {
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
        let fill = boolean(Op::Difference, &zone, &obst, FillRule::NonZero).unwrap();
        let tol = ArcTol::new(1000, Side::Inside);
        let t = Instant::now();
        let a = offset(&fill, -100_000, Join::Round, tol).unwrap();
        println!(
            "shrink: {:?} -> {} verts",
            t.elapsed(),
            a.iter().map(|p| p.vertex_count()).sum::<usize>()
        );
        let t = Instant::now();
        let b = offset(&a, 100_000, Join::Round, tol).unwrap();
        println!(
            "grow: {:?} -> {} verts",
            t.elapsed(),
            b.iter().map(|p| p.vertex_count()).sum::<usize>()
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
        let mut d = Vec::new();
        for _ in 0..20 {
            d = std::hint::black_box(
                boolean(Op::Difference, &zone, &obst, FillRule::NonZero).unwrap(),
            );
        }
        println!(
            "zone x20: {:?} (checksum {:x})",
            t.elapsed() / 20,
            checksum(&d)
        );
        return;
    }
    if std::env::args().any(|a| a == "incremental") {
        incremental_bench();
        return;
    }
    if std::env::args().any(|a| a == "bigdist") {
        let a = circle(0, 0, 50_000_000.0, 50_000);
        let b = circle(120_000_000, 0, 50_000_000.0, 50_000);
        let t = Instant::now();
        let c = distance(&a, &b).unwrap();
        println!(
            "distance between two 50k-vertex circles: {:?} -> {}",
            t.elapsed(),
            c.sq.distance_f64()
        );
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

/// Incremental zone refill: 100 mm zone − 5000 obstacles, then single and batched edits.
fn incremental_bench() {
    let mut s = 42u64;
    let zone = Ring::from([
        (0, 0),
        (100_000_000, 0),
        (100_000_000, 100_000_000),
        (0, 100_000_000),
    ]);
    let mut obst: Vec<Ring> = (0..5000)
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
    let mut z = ZoneFill::new(&zone, FillRule::NonZero).unwrap();
    for (i, o) in obst.iter().enumerate() {
        z.insert(i as u64, o).unwrap();
    }
    z.commit();
    println!("initial build: {:?}", t.elapsed());
    let t = Instant::now();
    let n = z.result().nodes.len();
    println!("first materialization: {:?} ({n} rings)", t.elapsed());
    let check = std::env::args().any(|a| a == "check");
    let verify = |z: &mut ZoneFill, obst: &[Ring]| {
        if check {
            let want = boolean(Op::Difference, &zone, obst, FillRule::NonZero).unwrap();
            assert_eq!(z.fill(), want);
        }
    };
    verify(&mut z, &obst);
    let rounds = 200;
    let rnd = |s: &mut u64| (lcg(s) % 100_000_000) as i64;
    // Single moves.
    let (mut tc, mut tf) = (std::time::Duration::ZERO, std::time::Duration::ZERO);
    for _ in 0..rounds {
        let id = (lcg(&mut s) % 5000) as usize;
        let r = circle(rnd(&mut s), rnd(&mut s), 300_000.0, 32);
        let t = Instant::now();
        z.update(id as u64, &r).unwrap();
        z.commit();
        tc += t.elapsed();
        let t = Instant::now();
        std::hint::black_box(z.result());
        tf += t.elapsed();
        obst[id] = r;
    }
    println!(
        "single move: commit {:?}, + result() {:?}, + fill() n/a",
        tc / rounds,
        tf / rounds
    );
    verify(&mut z, &obst);
    // Single insert + remove.
    let (mut ti, mut tr) = (std::time::Duration::ZERO, std::time::Duration::ZERO);
    for k in 0..rounds {
        let r = circle(rnd(&mut s), rnd(&mut s), 300_000.0, 32);
        let id = 1_000_000 + k as u64;
        let t = Instant::now();
        z.insert(id, &r).unwrap();
        z.commit();
        ti += t.elapsed();
        let t = Instant::now();
        z.remove(id);
        z.commit();
        tr += t.elapsed();
    }
    println!(
        "single insert: {:?}, single remove: {:?}",
        ti / rounds,
        tr / rounds
    );
    verify(&mut z, &obst);
    let t = Instant::now();
    for _ in 0..20 {
        std::hint::black_box(z.fill());
    }
    println!(
        "fill() materialization (cached tree -> PolygonSet): {:?}",
        t.elapsed() / 20
    );
    if std::env::args().any(|a| a == "noauto") {
        z.set_auto_rebuild(false);
    }
    for batch in [10usize, 100, 1000, 2500] {
        let reps = (2000 / batch).clamp(2, 50) as u32;
        let t = Instant::now();
        for _ in 0..reps {
            for _ in 0..batch {
                let id = (lcg(&mut s) % 5000) as usize;
                let r = circle(rnd(&mut s), rnd(&mut s), 300_000.0, 32);
                z.update(id as u64, &r).unwrap();
                obst[id] = r;
            }
            z.commit();
        }
        println!(
            "batch of {batch} moves: {:?} per batch (rebuilds so far: {})",
            t.elapsed() / reps,
            z.rebuild_count()
        );
    }
    verify(&mut z, &obst);
    let t = Instant::now();
    let d = boolean(Op::Difference, &zone, &obst, FillRule::NonZero).unwrap();
    println!("from-scratch boolean: {:?}", t.elapsed());
    std::hint::black_box(d);
}
