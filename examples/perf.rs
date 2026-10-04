use polyclip::*;

#[path = "../tests/common/corpus.rs"]
mod corpus;
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

fn best<T>(f: impl Fn() -> T) -> (std::time::Duration, T) {
    let mut d = std::time::Duration::MAX;
    let mut out = None;
    for _ in 0..3 {
        let t = Instant::now();
        out = Some(f());
        d = d.min(t.elapsed());
    }
    (d, out.unwrap())
}

/// POLYGON_LIB.md section 8: cadlab's opening workload on the real corpus.
fn cadlab() {
    let fill = corpus::load("cadlab_gnd_in1.pclp");
    let tol = ArcTol::new(5_000, Side::Inside);
    let w = 100_000;
    let ms = |d: std::time::Duration| d.as_secs_f64() * 1e3;
    let (d, shrunk) = best(|| offset(&fill, -w, Join::Round, tol).unwrap());
    println!("offset(-100 um): {:7.1} ms", ms(d));
    let (d, _) = best(|| union_all(&fill, FillRule::NonZero).unwrap());
    println!("  normalize input: {:7.1} ms", ms(d));
    let (d, grown) = best(|| offset(&shrunk, w, Join::Round, tol).unwrap());
    println!("offset(+100 um): {:7.1} ms", ms(d));
    let (d, opened) = best(|| opening(&fill, w, tol).unwrap());
    println!(
        "opening:         {:7.1} ms (equal: {})",
        ms(d),
        opened == grown
    );
    set_always_monolithic(true);
    let (d, one) = best(|| opening(&fill, w, tol).unwrap());
    set_always_monolithic(false);
    println!(
        "  in one piece:  {:7.1} ms (equal: {})",
        ms(d),
        opened == one
    );
    let spokes: Vec<Ring> = (0..40)
        .map(|k| {
            let (x, y) = (10_000_000 + k * 3_000_000, 50_000_000);
            Ring::from([
                (x, y),
                (x + 250_000, y),
                (x + 250_000, y + 1_000_000),
                (x, y + 1_000_000),
            ])
        })
        .collect();
    let spoked = |mono: bool| {
        Boolean::new()
            .subject(&opened, FillRule::NonZero)
            .subject(&spokes, FillRule::NonZero)
            .monolithic(mono)
            .execute()
            .unwrap()
    };
    let (d, u) = best(|| spoked(false));
    let (dm, um) = best(|| spoked(true));
    println!(
        "union + 40 spokes: {:7.1} ms (in one piece: {:.1} ms, equal: {})",
        ms(d),
        ms(dm),
        u == um
    );
}

/// Tracks: 45/90 degree polylines of a few segments, deterministic.
fn corpus_tracks(b: Rect, n: usize, seed: u64) -> Vec<Path> {
    let mut s = seed;
    let (w, h) = (b.width() as u64, b.height() as u64);
    (0..n)
        .map(|_| {
            let mut p = Point::new(
                b.min.x + (lcg(&mut s) % w) as i64,
                b.min.y + (lcg(&mut s) % h) as i64,
            );
            let mut pts = vec![p];
            for _ in 0..1 + lcg(&mut s) % 4 {
                let len = 500_000 + (lcg(&mut s) % 5_000_000) as i64;
                let (dx, dy) = [
                    (1, 0),
                    (1, 1),
                    (0, 1),
                    (-1, 1),
                    (-1, 0),
                    (-1, -1),
                    (0, -1),
                    (1, -1),
                ][(lcg(&mut s) % 8) as usize];
                p = Point::new(p.x + dx * len, p.y + dy * len);
                pts.push(p);
            }
            Path(pts)
        })
        .collect()
}

/// Pads: small axis-aligned rectangles, deterministic.
fn corpus_pads(b: Rect, n: usize, seed: u64) -> Vec<Ring> {
    let mut s = seed;
    let (w, h) = (b.width() as u64, b.height() as u64);
    (0..n)
        .map(|_| {
            let (x, y) = (
                b.min.x + (lcg(&mut s) % w) as i64,
                b.min.y + (lcg(&mut s) % h) as i64,
            );
            let (pw, ph) = (
                200_000 + (lcg(&mut s) % 1_500_000) as i64,
                200_000 + (lcg(&mut s) % 1_500_000) as i64,
            );
            Ring::from([(x, y), (x + pw, y), (x + pw, y + ph), (x, y + ph)])
        })
        .collect()
}

/// Operations other than booleans and offsets on the real corpus (DRC queries, path
/// clipping, fracture, triangulation, simplification, validation).
fn cadlab_ops() {
    let fill = corpus::load("cadlab_gnd_in1.pclp");
    let only: Vec<String> = std::env::args()
        .skip_while(|a| a != "cadlab-ops")
        .skip(1)
        .collect();
    let want = |k: &str| only.is_empty() || only.iter().any(|o| o == k);
    let us = |d: std::time::Duration| d.as_secs_f64() * 1e6;
    let ms = |d: std::time::Duration| d.as_secs_f64() * 1e3;
    let b = fill.bbox().unwrap();
    println!(
        "corpus: {} polygons ({} outer vertices), {} rings, {} vertices",
        fill.len(),
        fill[0].outer.len(),
        fill.iter().map(|p| 1 + p.holes.len()).sum::<usize>(),
        fill.iter().map(|p| p.vertex_count()).sum::<usize>()
    );
    let tracks = corpus_tracks(b, 500, 11);
    let pads = corpus_pads(b, 2000, 12);
    if want("clip1") {
        for _ in 0..300 {
            std::hint::black_box(clip_paths(&tracks[..1], &fill, FillRule::NonZero).unwrap());
        }
    }
    if want("clip") {
        let (d, r) = best(|| clip_paths(&tracks, &fill, FillRule::NonZero).unwrap());
        println!(
            "clip_paths 500 tracks: {:8.2} ms ({} in, {} out)",
            ms(d),
            r.inside.len(),
            r.outside.len()
        );
        set_always_monolithic(true);
        let (dm, rm) = best(|| clip_paths(&tracks, &fill, FillRule::NonZero).unwrap());
        set_always_monolithic(false);
        println!(
            "  in one piece:         {:8.2} ms (equal: {})",
            ms(dm),
            r == rm
        );
        let (d, _) = best(|| clip_paths(&tracks[..1], &fill, FillRule::NonZero).unwrap());
        println!("clip_paths 1 track:    {:8.2} ms", ms(d));
        let (d, _) = best(|| {
            tracks[..100]
                .iter()
                .map(|t| {
                    clip_paths(t, &fill, FillRule::NonZero)
                        .unwrap()
                        .inside
                        .len()
                })
                .sum::<usize>()
        });
        println!("clip_paths 100 x 1 track: {:8.2} ms/track", ms(d) / 100.0);
    }
    if want("drc") {
        let (d, n) = best(|| pads.iter().filter(|p| intersects(&fill, *p)).count());
        println!(
            "intersects:          {:8.2} us/pad ({n} hits)",
            us(d) / 2000.0
        );
        let (d, n) = best(|| {
            pads.iter()
                .filter(|p| distance_less_than(&fill, *p, 200_000))
                .count()
        });
        println!(
            "distance_less_than:  {:8.2} us/pad ({n} hits)",
            us(d) / 2000.0
        );
        let (d, n) = best(|| pads[..200].iter().filter(|p| contains(&fill, *p)).count());
        println!(
            "contains:            {:8.2} us/pad ({n} hits)",
            us(d) / 200.0
        );
        let (d, _) = best(|| {
            pads[..50]
                .iter()
                .map(|p| distance(&fill, p).unwrap().sq.to_f64())
                .sum::<f64>()
        });
        println!("distance:            {:8.2} us/pad", us(d) / 50.0);
        let (d, n) = best(|| {
            pads.iter()
                .filter(|p| locate(&fill, p[0]) == Location::Inside)
                .count()
        });
        println!(
            "locate:              {:8.2} us/pt ({n} inside)",
            us(d) / 2000.0
        );
        let (d, _) = best(|| area2(&fill));
        println!("area2:               {:8.2} ms", ms(d));
    }
    if want("prepared") || want("drc") {
        let t = Instant::now();
        let pz = Prepared::new(&fill);
        println!("Prepared::new:       {:8.2} ms", ms(t.elapsed()));
        // Pads centred in the holes (clearance checks) and the random ones.
        let hole_pads: Vec<Ring> = fill[0]
            .holes
            .iter()
            .take(2000)
            .map(|h| {
                let c = h.bbox().unwrap();
                let (x, y) = ((c.min.x + c.max.x) / 2, (c.min.y + c.max.y) / 2);
                Ring::from([
                    (x - 100_000, y - 100_000),
                    (x + 100_000, y - 100_000),
                    (x + 100_000, y + 100_000),
                    (x - 100_000, y + 100_000),
                ])
            })
            .collect();
        for (name, set) in [("random pads", &pads), ("pads in holes", &hole_pads)] {
            let n = set.len() as f64;
            println!("{name}:");
            let (d, a) = best(|| set.iter().filter(|p| pz.intersects(*p)).count());
            let (d0, b) = best(|| set[..100].iter().filter(|p| intersects(&fill, *p)).count());
            println!(
                "  intersects:         {:8.2} us/pad (free fn {:8.2} us) {a} {b}",
                us(d) / n,
                us(d0) / 100.0
            );
            let (d, a) = best(|| {
                set.iter()
                    .filter(|p| pz.distance_less_than(*p, 200_000))
                    .count()
            });
            let (d0, b) = best(|| {
                set[..100]
                    .iter()
                    .filter(|p| distance_less_than(&fill, *p, 200_000))
                    .count()
            });
            println!(
                "  distance_less_than: {:8.2} us/pad (free fn {:8.2} us) {a} {b}",
                us(d) / n,
                us(d0) / 100.0
            );
            let (d, a) = best(|| set.iter().filter(|p| pz.contains(*p)).count());
            let (d0, b) = best(|| set[..20].iter().filter(|p| contains(&fill, *p)).count());
            println!(
                "  contains:           {:8.2} us/pad (free fn {:8.2} us) {a} {b}",
                us(d) / n,
                us(d0) / 20.0
            );
            let (d, _) = best(|| {
                set.iter()
                    .map(|p| pz.distance(p).unwrap().sq.to_f64())
                    .sum::<f64>()
            });
            let (d0, _) = best(|| {
                set[..20]
                    .iter()
                    .map(|p| distance(&fill, p).unwrap().sq.to_f64())
                    .sum::<f64>()
            });
            println!(
                "  distance:           {:8.2} us/pad (free fn {:8.2} us)",
                us(d) / n,
                us(d0) / 20.0
            );
            let (d, a) = best(|| {
                set.iter()
                    .filter(|p| pz.locate(p[0]) == Location::Inside)
                    .count()
            });
            println!("  locate:             {:8.2} us/pt {a}", us(d) / n);
        }
    }
    if want("fracture") {
        let (d, r) = best(|| fracture_set(&fill).unwrap());
        println!(
            "fracture:            {:8.2} ms ({} verts)",
            ms(d),
            r.iter().map(|r| r.len()).sum::<usize>()
        );
    }
    if want("triangulate") {
        let (d, t) = best(|| triangulate_set(&fill).unwrap());
        let hash = |t: &Triangulation| {
            use std::hash::{Hash, Hasher};
            let mut h = std::collections::hash_map::DefaultHasher::new();
            t.hash(&mut h);
            h.finish()
        };
        println!(
            "triangulate:         {:8.2} ms ({} triangles, checksum {:x})",
            ms(d),
            t.triangles.len(),
            hash(&t)
        );
        let (d, t) = best(|| triangulate_delaunay(&fill[0]).unwrap());
        println!(
            "triangulate_delaunay:{:8.2} ms (checksum {:x})",
            ms(d),
            hash(&t)
        );
    }
    if want("simplify") {
        let (d, r) = best(|| simplify_polygons(&fill, 5_000));
        println!(
            "simplify(5 um):      {:8.2} ms ({} verts)",
            ms(d),
            r.iter().map(|p| p.vertex_count()).sum::<usize>()
        );
    }
    if want("validate") {
        let (d, _) = best(|| validate_set(&fill).unwrap());
        println!("validate_set:        {:8.2} ms", ms(d));
        let (d, _) = best(|| check_canonical(&fill, true).unwrap());
        println!("check_canonical:     {:8.2} ms", ms(d));
    }
    if want("zonefill") {
        // The pour as the zone, pads as obstacles.
        let t = Instant::now();
        let mut z = ZoneFill::new(&fill, FillRule::NonZero).unwrap();
        for (i, p) in pads[..500].iter().enumerate() {
            z.insert(i as u64, p).unwrap();
        }
        z.commit();
        std::hint::black_box(z.result());
        println!("ZoneFill build:      {:8.2} ms", ms(t.elapsed()));
        let mut s = 99u64;
        let t = Instant::now();
        for _ in 0..100 {
            let id = lcg(&mut s) % 500;
            let p = &pads[500 + (lcg(&mut s) % 1500) as usize];
            z.update(id, p).unwrap();
            z.commit();
            std::hint::black_box(z.result());
        }
        println!("ZoneFill move+result:{:8.2} us", us(t.elapsed()) / 100.0);
    }
}

fn main() {
    if std::env::args().any(|a| a == "cadlab-ops") {
        cadlab_ops();
        return;
    }
    if std::env::args().any(|a| a == "cadlab") {
        cadlab();
        return;
    }
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
