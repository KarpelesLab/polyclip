//! Round-join offsets of pour-shaped regions: shrink, grow and the opening (shrink, then
//! grow), as cadlab's zone fills run them for the minimum copper width.
//!
//! ```text
//! cargo run --release --example offset_bench                 # synthetic pours
//! cargo run --release --example offset_bench -- FILE [REPS]  # recorded openings
//! cargo run --release --example offset_bench -- synthetic 1  # synthetic, timed once
//! cargo run --release --example offset_bench -- FILE 1 11    # call 11 of the file alone
//! ```
//!
//! Without a file the inputs are generated (fixed seed): two layers of a board of 100 mm
//! with a pour over the whole board, from which the pads of footprints, vias and tracks
//! grown by their clearance are taken, each poured four times while tracks are added. That
//! gives the statistics of recorded pours (hundreds of rings, 10^4 to 10^5 vertices, most
//! of them on arcs of 5 um tolerance, the same pad many times over, most of a pour as it
//! was at the fill before; coordinates in nanometres).
//!
//! A file holds recorded calls, as little-endian `i64`: for each call the offset `d`, the
//! arc tolerance and the number of polygons, then for each polygon its number of rings
//! (outer first), and for each ring its number of points and their `x`, `y`.
//!
//! With `OFFSET_BENCH_MEMO=<MiB>` the openings run with a [`Memo`] of that size kept over
//! the calls, and with `OFFSET_BENCH_VERIFY` every reused result is computed and compared.
//!
//! Every call is timed `REPS` times (default 3, the best kept). The last line has a hash of
//! every vertex of every result in order: two builds that print the same hash returned the
//! same polygons.

use polyclip::*;
use std::time::{Duration, Instant};

struct Case {
    d: i64,
    tol: i64,
    input: PolygonSet,
}

fn load(path: &str) -> Vec<Case> {
    let bytes = std::fs::read(path).expect("cannot read the dump");
    let mut words = bytes
        .as_chunks::<8>()
        .0
        .iter()
        .map(|c| i64::from_le_bytes(*c));
    let mut next = || words.next().expect("truncated dump");
    let mut cases = Vec::new();
    let total = bytes.len() / 8;
    let mut used = 0usize;
    while used < total {
        let (d, tol, np) = (next(), next(), next());
        used += 3;
        let mut input: Vec<Polygon> = Vec::with_capacity(np as usize);
        for _ in 0..np {
            let nr = next();
            used += 1;
            let mut rings: Vec<Ring> = Vec::with_capacity(nr as usize);
            for _ in 0..nr {
                let n = next();
                used += 1 + 2 * n as usize;
                rings.push((0..n).map(|_| Point::new(next(), next())).collect());
            }
            let mut rings = rings.into_iter();
            let outer = rings.next().expect("polygon without rings");
            input.push(Polygon::new(outer, rings.collect()));
        }
        cases.push(Case { d, tol, input });
    }
    cases
}

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

const MM: i64 = 1_000_000;
/// Side of the synthetic board.
const SIDE: i64 = 100 * MM;

/// The copper of a layer other than its pour: pads of footprints (rows of equal pads at a
/// pitch, as integrated circuits and connectors have them, pairs of pads of passives), round
/// vias, and tracks of a few segments along the eight directions.
struct Layer {
    pads: Vec<Ring>,
    tracks: Vec<Path>,
}

fn layer(seed: u64, footprints: usize, vias: usize, tracks: usize) -> Layer {
    let mut s = Lcg(seed);
    let tol = ArcTol::new(5_000, Side::Outside);
    let rect = |x: i64, y: i64, w: i64, h: i64| {
        Ring::from([
            (x - w, y - h),
            (x + w, y - h),
            (x + w, y + h),
            (x - w, y + h),
        ])
    };
    let mut pads: Vec<Ring> = Vec::new();
    for _ in 0..footprints {
        let (x, y) = (
            s.below(SIDE - 20 * MM) + 5 * MM,
            s.below(SIDE - 20 * MM) + 5 * MM,
        );
        match s.below(3) {
            // Two rows of pads, the rows along x or along y.
            0 | 1 => {
                let n = 4 + s.below(21);
                let pitch = [500_000, 650_000, 1_270_000][s.below(3) as usize];
                let (w, h) = (pitch * 3 / 10, 400_000 + s.below(600_000));
                let rows = 2 * MM + s.below(6 * MM);
                let turned = s.below(2) == 0;
                for k in 0..n {
                    for row in [0, rows] {
                        pads.push(if turned {
                            rect(x + row, y + k * pitch, h, w)
                        } else {
                            rect(x + k * pitch, y + row, w, h)
                        });
                    }
                }
            }
            // A passive: two pads.
            _ => {
                let (w, h) = (300_000 + s.below(300_000), 400_000 + s.below(400_000));
                pads.push(rect(x, y, w, h));
                pads.push(rect(x + 3 * w, y, w, h));
            }
        }
    }
    // Vias of one size: the same polygon moved.
    let via = Circle::new(Point::new(0, 0), 300_000).to_ring(tol).unwrap();
    for _ in 0..vias {
        let (x, y) = (s.below(SIDE), s.below(SIDE));
        pads.push(via.0.iter().map(|p| Point::new(p.x + x, p.y + y)).collect());
    }
    let dirs: [(i64, i64); 8] = [
        (1, 0),
        (1, 1),
        (0, 1),
        (-1, 1),
        (-1, 0),
        (-1, -1),
        (0, -1),
        (1, -1),
    ];
    let mut paths: Vec<Path> = Vec::with_capacity(tracks);
    for _ in 0..tracks {
        let mut p = Point::new(s.below(SIDE), s.below(SIDE));
        let mut pts = vec![p];
        let mut dir = s.below(8) as usize;
        for _ in 0..2 + s.below(4) {
            let len = MM / 2 + s.below(6 * MM);
            p = Point::new(
                (p.x + dirs[dir].0 * len).clamp(0, SIDE),
                (p.y + dirs[dir].1 * len).clamp(0, SIDE),
            );
            pts.push(p);
            dir = (dir + if s.below(2) == 0 { 1 } else { 7 }) % 8;
        }
        paths.push(pts.into_iter().collect());
    }
    Layer {
        pads,
        tracks: paths,
    }
}

/// The pour of a layer with the first `tracks` of its tracks: the board less the pads and
/// those tracks (0.15 mm wide), each grown by a clearance of 0.2 mm with arcs of 5 um.
fn pour(l: &Layer, tracks: usize) -> PolygonSet {
    let tol = ArcTol::new(5_000, Side::Outside);
    let mut obstacles = offset_paths(
        &l.tracks[..tracks.min(l.tracks.len())],
        75_000 + 200_000,
        Join::Round,
        EndCap::Round,
        tol,
    )
    .unwrap();
    obstacles.extend(offset(&l.pads, 200_000, Join::Round, tol).unwrap());
    let board = Ring::from([(0, 0), (SIDE, 0), (SIDE, SIDE), (0, SIDE)]);
    boolean(Op::Difference, &board, &obstacles, FillRule::NonZero).unwrap()
}

/// Two layers (pads and a few tracks; many tracks), each poured four times while tracks are
/// added to it, as a router refills its zones.
fn synthetic() -> Vec<Case> {
    let mut cases = Vec::new();
    // (seed, footprints, vias, tracks at first, tracks added by fill, d)
    for (seed, footprints, vias, tracks, step, d) in [
        (1, 60, 250, 150, 25, 100_000),
        (2, 30, 500, 900, 60, 127_000),
    ] {
        let l = layer(seed, footprints, vias, tracks + 3 * step);
        for fill in 0..4 {
            cases.push(Case {
                d,
                tol: 5_000,
                input: pour(&l, tracks + fill * step),
            });
        }
    }
    cases
}

fn mix(h: &mut u64, x: i64) {
    for b in x.to_le_bytes() {
        *h ^= b as u64;
        *h = h.wrapping_mul(0x100000001b3);
    }
}

fn hash(h: &mut u64, s: &PolygonSet) {
    mix(h, -3);
    for p in s {
        mix(h, -1);
        for r in p.rings() {
            mix(h, -2);
            for q in &r.0 {
                mix(h, q.x);
                mix(h, q.y);
            }
        }
    }
}

fn best<T>(reps: usize, f: impl Fn() -> T) -> (Duration, T) {
    let mut d = Duration::MAX;
    let mut out = None;
    for _ in 0..reps.max(1) {
        let t = Instant::now();
        out = Some(f());
        d = d.min(t.elapsed());
    }
    (d, out.unwrap())
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let cases = match args.first() {
        Some(path) if path != "synthetic" => load(path),
        _ => synthetic(),
    };
    let reps: usize = args.get(1).and_then(|r| r.parse().ok()).unwrap_or(3);
    // A third argument times that call alone.
    let only: Option<usize> = args.get(2).and_then(|r| r.parse().ok());
    let verbose = std::env::var_os("OFFSET_BENCH_VERBOSE").is_some()
        || args.first().is_none_or(|a| a == "synthetic");
    let ms = |d: Duration| d.as_secs_f64() * 1e3;
    // OFFSET_BENCH_MEMO=<MiB>: the openings with a memo of that size kept over the calls, as
    // a caller refilling the same zones has it. OFFSET_BENCH_VERIFY: every result that is
    // reused is computed as well and compared.
    let memo: Option<Memo> = std::env::var("OFFSET_BENCH_MEMO")
        .ok()
        .and_then(|m| m.parse::<usize>().ok())
        .map(|mib| Memo::new(mib << 20));
    set_memo_verify(std::env::var_os("OFFSET_BENCH_VERIFY").is_some());
    let mut h: u64 = 0xcbf29ce484222325;
    let (mut t_in, mut t_out, mut t_open) = (Duration::ZERO, Duration::ZERO, Duration::ZERO);
    let (mut verts, mut rings) = (0usize, 0usize);
    let mut largest = (0usize, 0usize, Duration::ZERO);
    for (i, c) in cases.iter().enumerate() {
        if only.is_some_and(|o| o != i) {
            continue;
        }
        let tol = ArcTol::new(c.tol, Side::Inside);
        let nv: usize = c.input.iter().map(|p| p.vertex_count()).sum();
        let nr: usize = c.input.iter().map(|p| 1 + p.holes.len()).sum();
        let (d_in, shrunk) = best(reps, || offset(&c.input, -c.d, Join::Round, tol).unwrap());
        let (d_out, grown) = best(reps, || offset(&shrunk, c.d, Join::Round, tol).unwrap());
        let (d_open, opened) = match &memo {
            // Timed once: a second run would find everything of the first.
            Some(m) => best(1, || m.opening(&c.input, c.d, tol).unwrap()),
            None => best(reps, || opening(&c.input, c.d, tol).unwrap()),
        };
        assert!(opened == grown, "the opening is not shrink, then grow");
        hash(&mut h, &shrunk);
        hash(&mut h, &opened);
        t_in += d_in;
        t_out += d_out;
        t_open += d_open;
        verts += nv;
        rings += nr;
        if nv > largest.0 {
            largest = (nv, nr, d_open);
        }
        if verbose {
            println!(
                "call {i}: d {} um, {nv} vertices in {nr} rings: shrink {:.2} ms, grow {:.2} ms, opening {:.2} ms -> {} vertices",
                c.d / 1000,
                ms(d_in),
                ms(d_out),
                ms(d_open),
                opened.iter().map(|p| p.vertex_count()).sum::<usize>()
            );
        }
    }
    println!(
        "{} calls, {verts} vertices in {rings} rings: shrink {:.1} ms, grow {:.1} ms, opening {:.1} ms ({:.2} us a vertex); largest {} vertices in {} rings {:.1} ms; hash {h:016x}",
        cases.len(),
        ms(t_in),
        ms(t_out),
        ms(t_open),
        t_open.as_secs_f64() * 1e6 / verts.max(1) as f64,
        largest.0,
        largest.1,
        ms(largest.2),
    );
    if let Some(m) = &memo {
        let s = m.stats();
        println!(
            "memo: {} clusters, {:.1} MiB, {} found, {} computed",
            s.entries,
            s.bytes as f64 / (1 << 20) as f64,
            s.hits,
            s.misses
        );
    }
}
