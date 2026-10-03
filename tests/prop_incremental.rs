//! The incremental zone fill must equal the from-scratch boolean after every edit.
//!
//! Random sequences of insert / update / remove / clear / batches on a [`ZoneFill`], with
//! obstacles that overlap, touch, coincide, cross the zone boundary, lie far away, are
//! degenerate, carry tags, or are thin diagonals that force heavy snap rounding. After every
//! step the result is compared (`==` on the canonical tree, and on tagged polygons) with
//! `Boolean` run from scratch on the current obstacle set.

use polyclip::*;
use proptest::prelude::*;
use std::collections::BTreeMap;

fn rect(x: i64, y: i64, w: i64, h: i64) -> Vec<Point> {
    vec![
        Point::new(x, y),
        Point::new(x + w, y),
        Point::new(x + w, y + h),
        Point::new(x, y + h),
    ]
}

/// An obstacle: one or more rings with per-edge tags.
#[derive(Clone, Debug)]
struct Obst(Vec<TaggedRing>);

fn tagged(pts: Vec<Point>, tag: u64, per_edge: bool) -> TaggedRing {
    let n = pts.len();
    TaggedRing {
        points: pts,
        tags: (0..n as u64)
            .map(|i| if per_edge { tag * 16 + i % 3 } else { tag })
            .collect(),
    }
}

/// Obstacle geometry at scale `range` (the zone is roughly `[0, range]^2`).
fn obstacle(range: i64) -> impl Strategy<Value = Obst> {
    let r = range;
    let random_poly = (
        prop::collection::vec((-r / 4..r + r / 4, -r / 4..r + r / 4), 0..8),
        0u64..4,
        any::<bool>(),
    )
        .prop_map(|(v, t, pe)| {
            Obst(vec![tagged(
                v.into_iter().map(Point::from).collect(),
                t,
                pe,
            )])
        });
    // Axis-aligned rectangles on a coarse grid: touching, sharing edges, coinciding.
    let grid_rect = (0..8i64, 0..8i64, 1..4i64, 1..4i64, -1..9i64, 0u64..3).prop_map(
        move |(x, y, w, h, shift, t)| {
            let s = (r / 8).max(1);
            let off = if shift < 0 { -s / 2 } else { 0 };
            Obst(vec![tagged(
                rect(x * s + off, y * s, w * s, h * s),
                t,
                false,
            )])
        },
    );
    // Thin diagonal slivers (crossings far from integer points: heavy rounding).
    let diag = (0..r, 0..r, 1..r.max(2), 1i64..4, 0u64..3).prop_map(move |(x, y, len, w, t)| {
        Obst(vec![tagged(
            vec![
                Point::new(x, y),
                Point::new(x + len, y + len * 2 / 3),
                Point::new(x + len - w, y + len * 2 / 3 + w),
                Point::new(x - w, y + w),
            ],
            t,
            true,
        )])
    });
    // A polygon with a hole (two rings), possibly crossing the zone boundary.
    let holed = (-r / 4..r, -r / 4..r, 4..r.max(5), 0u64..3).prop_map(move |(x, y, s, t)| {
        let mut hole = rect(x + s / 4, y + s / 4, s / 2, s / 2);
        hole.reverse();
        Obst(vec![
            tagged(rect(x, y, s, s), t, false),
            tagged(hole, t + 7, false),
        ])
    });
    // Tiny triangles (a few units).
    let tiny = (0..r, 0..r, 0u64..3).prop_map(|(x, y, t)| {
        Obst(vec![tagged(
            vec![
                Point::new(x, y),
                Point::new(x + 3, y + 1),
                Point::new(x + 1, y + 2),
            ],
            t,
            false,
        )])
    });
    // Far away or huge shapes (other scales; long segments and sparse grid cells).
    let far = (0..3i64, 0u64..3).prop_map(move |(k, t)| {
        let big = [r * 1000, 1 << 39, -(1 << 39)][k as usize];
        Obst(vec![tagged(rect(big, big / 2, r, r * 3), t, false)])
    });
    let spanning = (0..r, 0u64..3).prop_map(move |(y, t)| {
        // A long sliver across the whole zone and far beyond.
        Obst(vec![tagged(
            vec![
                Point::new(-(1 << 38), y),
                Point::new(1 << 38, y + 1),
                Point::new(1 << 38, y + 3),
            ],
            t,
            false,
        )])
    });
    prop_oneof![
        5 => random_poly,
        4 => grid_rect,
        3 => diag,
        2 => holed,
        2 => tiny,
        1 => far,
        1 => spanning,
    ]
}

#[derive(Clone, Debug)]
enum Step {
    Insert(u64, Obst),
    Update(u64, Obst),
    Remove(u64),
    /// Re-insert a copy of another obstacle's geometry (duplicates coincide exactly).
    Duplicate(u64, u64),
    Clear,
    /// Several edits before the next comparison (one batch).
    Batch(Vec<(u64, Option<Obst>)>),
}

fn step(range: i64) -> impl Strategy<Value = Step> {
    let id = 0u64..10;
    prop_oneof![
        6 => (id.clone(), obstacle(range)).prop_map(|(i, o)| Step::Insert(i, o)),
        3 => (id.clone(), obstacle(range)).prop_map(|(i, o)| Step::Update(i, o)),
        3 => id.clone().prop_map(Step::Remove),
        1 => (id.clone(), id.clone()).prop_map(|(a, b)| Step::Duplicate(a, b)),
        1 => Just(Step::Clear),
        2 => prop::collection::vec((id, prop::option::of(obstacle(range))), 1..5)
            .prop_map(Step::Batch),
    ]
}

fn zone(range: i64) -> impl Strategy<Value = Vec<TaggedRing>> {
    let r = range;
    prop_oneof![
        3 => Just(vec![tagged(rect(0, 0, r, r), 100, false)]),
        // Square with a square hole.
        1 => Just({
            let mut hole = rect(r / 3, r / 3, r / 3, r / 3);
            hole.reverse();
            vec![tagged(rect(0, 0, r, r), 100, false), tagged(hole, 101, false)]
        }),
        // Arbitrary (self-intersecting) zone outline.
        1 => prop::collection::vec((0..r, 0..r), 3..9).prop_map(|v| {
            vec![tagged(v.into_iter().map(Point::from).collect(), 102, true)]
        }),
        // Empty zone.
        1 => Just(Vec::new()),
    ]
}

fn rule() -> impl Strategy<Value = FillRule> {
    prop_oneof![
        4 => Just(FillRule::NonZero),
        1 => Just(FillRule::EvenOdd),
        1 => Just(FillRule::Positive),
        1 => Just(FillRule::Negative)
    ]
}

fn reference(zone: &[TaggedRing], obst: &BTreeMap<u64, Obst>, rule: FillRule) -> PolyTree {
    let clip: Vec<TaggedRing> = obst.values().flat_map(|o| o.0.iter().cloned()).collect();
    Boolean::new()
        .subject(zone, rule)
        .clip(&clip, rule)
        .op(Op::Difference)
        .execute_tree()
        .unwrap()
}

fn run(zone: Vec<TaggedRing>, rule: FillRule, auto: bool, steps: Vec<Step>) {
    let mut z = ZoneFill::new(&zone, rule).unwrap();
    z.set_auto_rebuild(auto);
    let mut cur: BTreeMap<u64, Obst> = BTreeMap::new();
    for (k, s) in steps.into_iter().enumerate() {
        match s {
            Step::Insert(i, o) => {
                assert_eq!(z.insert(i, &o.0).unwrap(), cur.contains_key(&i));
                cur.insert(i, o);
            }
            Step::Update(i, o) => {
                assert_eq!(z.update(i, &o.0).unwrap(), cur.contains_key(&i));
                if cur.contains_key(&i) {
                    cur.insert(i, o);
                }
            }
            Step::Remove(i) => {
                assert_eq!(z.remove(i), cur.remove(&i).is_some());
            }
            Step::Duplicate(a, b) => {
                if let Some(o) = cur.get(&b).cloned() {
                    z.insert(a, &o.0).unwrap();
                    cur.insert(a, o);
                }
            }
            Step::Clear => {
                z.clear();
                cur.clear();
            }
            Step::Batch(v) => {
                for (i, o) in v {
                    match o {
                        Some(o) => {
                            z.insert(i, &o.0).unwrap();
                            cur.insert(i, o);
                        }
                        None => {
                            z.remove(i);
                            cur.remove(&i);
                        }
                    }
                }
            }
        }
        assert_eq!(z.len(), cur.len());
        assert_eq!(z.ids(), cur.keys().copied().collect::<Vec<_>>());
        let want = reference(&zone, &cur, rule);
        if k % 2 == 0 {
            // Flattened results computed directly (before the tree is cached).
            assert_eq!(z.fill(), want.to_polygon_set(), "set, step {k}");
            assert_eq!(
                z.fill_tagged(),
                want.to_tagged_polygons(),
                "tagged, step {k}"
            );
        }
        assert_eq!(z.fill_tree(), want, "step {k}");
        // The internal state of every stage matches a from-scratch recomputation too.
        if let Err(e) = z.verify() {
            panic!("step {k}: {e}");
        }
        if k % 2 == 1 {
            // Flattened from the cached tree.
            assert_eq!(
                z.fill_tagged(),
                want.to_tagged_polygons(),
                "tagged, step {k}"
            );
            assert_eq!(z.fill(), want.to_polygon_set(), "set, step {k}");
        }
    }
}

type Case = (Vec<TaggedRing>, FillRule, bool, Vec<Step>);

fn case() -> impl Strategy<Value = Case> {
    prop_oneof![Just(12i64), Just(100), Just(5_000), Just(1_000_000)].prop_flat_map(|r| {
        (
            zone(r),
            rule(),
            prop::bool::weighted(0.2),
            prop::collection::vec(step(r), 1..25),
        )
    })
}

proptest! {
    #[test]
    fn incremental_equals_from_scratch((zone, rule, auto, steps) in case()) {
        run(zone, rule, auto, steps);
    }
}

fn lcg(s: &mut u64) -> u64 {
    *s = s
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    *s >> 33
}

/// Deterministic random sequences on small grids, where snap rounding cascades are common
/// (many crossings within a few units).
#[test]
fn dense_rounding_sequences() {
    for seed in 0..300u64 {
        let mut s = seed * 7919 + 1;
        let range = [20i64, 100, 1000, 100_000][(seed % 4) as usize];
        let zone = Ring::from([(0, 0), (range, 0), (range, range), (0, range)]);
        let mut z = ZoneFill::new(&zone, FillRule::NonZero).unwrap();
        z.set_auto_rebuild(false);
        let mut cur: BTreeMap<u64, Ring> = BTreeMap::new();
        for step in 0..30 {
            let id = lcg(&mut s) % 8;
            if lcg(&mut s) % 3 < 2 {
                let n = 3 + (lcg(&mut s) % 5) as usize;
                let cx = (lcg(&mut s) % range as u64) as i64;
                let cy = (lcg(&mut s) % range as u64) as i64;
                let r = 1 + (lcg(&mut s) % (range as u64 / 3)) as i64;
                let ring: Ring = (0..n)
                    .map(|_| {
                        Point::new(
                            cx + (lcg(&mut s) % (2 * r as u64)) as i64 - r,
                            cy + (lcg(&mut s) % (2 * r as u64)) as i64 - r,
                        )
                    })
                    .collect();
                z.insert(id, &ring).unwrap();
                cur.insert(id, ring);
            } else {
                z.remove(id);
                cur.remove(&id);
            }
            let obst: Vec<Ring> = cur.values().cloned().collect();
            let want = Boolean::new()
                .subject(&zone, FillRule::NonZero)
                .clip(&obst, FillRule::NonZero)
                .op(Op::Difference)
                .execute_tree()
                .unwrap();
            assert_eq!(z.fill_tree(), want, "seed {seed} step {step}");
            if let Err(e) = z.verify() {
                panic!("seed {seed} step {step}: {e}");
            }
        }
        // Every update took the incremental path (no consistency fallback).
        assert_eq!(z.rebuild_count(), 1, "seed {seed}");
    }
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

/// A zone-fill workload in miniature: many round obstacles, then small batches of moves.
#[test]
fn zone_with_circles_moves() {
    let side = 10_000_000i64;
    let zone = Ring::from([(0, 0), (side, 0), (side, side), (0, side)]);
    let mut s = 7u64;
    let rnd = |s: &mut u64| (lcg(s) % side as u64) as i64;
    let mut obst: Vec<Ring> = (0..300)
        .map(|_| circle(rnd(&mut s), rnd(&mut s), 300_000.0, 24))
        .collect();
    let mut z = ZoneFill::new(&zone, FillRule::NonZero).unwrap();
    for (i, o) in obst.iter().enumerate() {
        z.insert(i as u64, o).unwrap();
    }
    assert_eq!(
        z.fill(),
        boolean(Op::Difference, &zone, &obst, FillRule::NonZero).unwrap()
    );
    let rebuilds = z.rebuild_count();
    for round in 0..60 {
        for _ in 0..1 + round % 4 {
            let i = (lcg(&mut s) % obst.len() as u64) as usize;
            let r = circle(rnd(&mut s), rnd(&mut s), 300_000.0, 24);
            z.update(i as u64, &r).unwrap();
            obst[i] = r;
        }
        assert_eq!(
            z.fill(),
            boolean(Op::Difference, &zone, &obst, FillRule::NonZero).unwrap(),
            "round {round}"
        );
    }
    assert_eq!(
        z.rebuild_count(),
        rebuilds,
        "small batches stay incremental"
    );
}

#[test]
fn send_sync() {
    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<ZoneFill>();
}
