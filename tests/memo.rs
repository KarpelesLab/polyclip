//! Results reused for clusters of rings that repeat, within one operation (the same rings at
//! another place) and across operations with a [`Memo`], must be the results computed in one
//! piece, for every op, fill rule and input.

use polyclip::*;
use proptest::prelude::*;

struct Lcg(u64);

impl Lcg {
    fn below(&mut self, n: i64) -> i64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 33) % n as u64) as i64
    }
}

/// A few ring shapes around the origin: simple, self-crossing, with collinear vertices.
fn shape(kind: i64, s: i64) -> Vec<(i64, i64)> {
    match kind {
        0 => vec![(0, 0), (4 * s, 0), (4 * s, 3 * s), (0, 3 * s)],
        1 => vec![(0, 0), (4 * s, 3 * s), (4 * s, 0), (0, 3 * s)],
        2 => vec![
            (0, 0),
            (2 * s, 0),
            (4 * s, 0),
            (4 * s, 4 * s),
            (3 * s, s),
            (s, 5 * s),
        ],
        3 => vec![
            (0, 0),
            (5 * s, s),
            (s, 2 * s),
            (5 * s, 3 * s),
            (0, 4 * s),
            (3 * s, 2 * s),
        ],
        _ => vec![
            (0, 0),
            (3 * s, 5 * s),
            (6 * s, 0),
            (0, 3 * s),
            (6 * s, 3 * s),
        ],
    }
}

fn at(pts: &[(i64, i64)], x: i64, y: i64, reverse: bool) -> Ring {
    let mut r: Vec<Point> = pts.iter().map(|&(a, b)| Point::new(a + x, b + y)).collect();
    if reverse {
        r.reverse();
    }
    r.into_iter().collect()
}

/// `n` rings of a few shapes on a grid of `pitch`, with pairs of rings and stray rings close
/// by, touching or overlapping: many clusters that are one another moved.
fn scene(seed: u64, n: usize, pitch: i64, s: i64) -> Vec<Ring> {
    let mut rnd = Lcg(seed);
    let side = (n as f64).sqrt().ceil() as i64;
    let mut rings = Vec::new();
    for i in 0..n as i64 {
        let (x, y) = (pitch * (i % side), pitch * (i / side));
        let kind = rnd.below(5);
        rings.push(at(&shape(kind, s), x, y, rnd.below(4) == 0));
        match rnd.below(6) {
            // A second ring at a set place next to the first.
            0 => rings.push(at(&shape((kind + 1) % 5, s), x + 6 * s + 2, y, false)),
            1 => rings.push(at(&shape(0, s), x + 2 * s, y + s, true)),
            // A stray ring anywhere near.
            2 => rings.push(at(
                &shape(rnd.below(5), s),
                x + rnd.below(8 * s) - 4 * s,
                y + rnd.below(8 * s) - 4 * s,
                false,
            )),
            _ => {}
        }
    }
    rings
}

fn fill(k: u8) -> FillRule {
    [
        FillRule::NonZero,
        FillRule::EvenOdd,
        FillRule::Positive,
        FillRule::Negative,
    ][k as usize % 4]
}

fn op(k: u8) -> Op {
    [Op::Union, Op::Intersection, Op::Difference, Op::Xor][k as usize % 4]
}

/// Runs `f` with every reused result computed as well and compared (a panic otherwise).
/// The switch is the process's and stays on: the tests of this file run side by side.
fn verified<T>(f: impl FnOnce() -> T) -> T {
    set_memo_verify(true);
    f()
}

#[test]
fn moved_clusters_identical() {
    for seed in 0..6u64 {
        let rings = scene(seed, 400, 97, 7);
        for k in 0..4 {
            let run = |mono: bool| {
                Boolean::new()
                    .subject(&rings, fill(k))
                    .monolithic(mono)
                    .execute_tree()
                    .unwrap()
            };
            assert!(verified(|| run(false)) == run(true), "seed {seed} rule {k}");
        }
    }
}

#[test]
fn memo_booleans_identical() {
    let memo = Memo::new(8 << 20);
    for seed in 0..4u64 {
        let subject = scene(seed, 300, 97, 7);
        let clip = scene(seed + 100, 300, 97, 5);
        for k in 0..16u8 {
            let run = |memo: Option<&Memo>| {
                let mut b = Boolean::new()
                    .subject(&subject, fill(k))
                    .clip(&clip, fill(k / 4 + 1))
                    .op(op(k))
                    .keep_collinear(k % 3 == 0);
                b.set_memo(memo);
                b.execute_tree().unwrap()
            };
            let plain = Boolean::new()
                .subject(&subject, fill(k))
                .clip(&clip, fill(k / 4 + 1))
                .op(op(k))
                .keep_collinear(k % 3 == 0)
                .monolithic(true)
                .execute_tree()
                .unwrap();
            let first = verified(|| run(Some(&memo)));
            let again = verified(|| run(Some(&memo)));
            assert!(first == plain, "seed {seed} case {k}: first run");
            assert!(again == plain, "seed {seed} case {k}: second run");
        }
    }
    let s = memo.stats();
    assert!(s.hits > 0 && s.misses > 0 && s.entries > 0, "{s:?}");
}

#[test]
fn memo_offsets_identical() {
    let memo = Memo::new(16 << 20);
    let tol = ArcTol::new(2, Side::Inside);
    for seed in 0..3u64 {
        // A scene, then the same with a part of it replaced: most clusters as they were.
        let a = scene(seed, 300, 130, 9);
        let mut b = a.clone();
        b.truncate(a.len() * 9 / 10);
        b.extend(scene(seed + 50, 30, 130, 9).into_iter().map(|r| {
            r.0.iter()
                .map(|p| Point::new(p.x + 40, p.y - 3000))
                .collect::<Ring>()
        }));
        for input in [&a, &b] {
            let plain = |d: i64, join: Join| {
                set_always_monolithic(true);
                let out = offset(input, d, join, tol).unwrap();
                set_always_monolithic(false);
                out
            };
            for (d, join) in [
                (6, Join::Round),
                (-4, Join::Round),
                (5, Join::Miter { limit: 2.0 }),
                (-3, Join::Bevel),
            ] {
                let got = verified(|| memo.offset(input, d, join, tol).unwrap());
                assert!(got == plain(d, join), "seed {seed} offset {d} {join:?}");
                let tree = verified(|| memo.offset_tagged(input, d, join, tol, 7).unwrap());
                assert!(
                    tree == offset_tagged(input, d, join, tol, 7).unwrap(),
                    "seed {seed} tagged offset {d} {join:?}"
                );
            }
            let opened = verified(|| memo.opening(input, 5, tol).unwrap());
            assert!(
                opened == opening(input, 5, tol).unwrap(),
                "seed {seed} opening"
            );
            let closed = verified(|| memo.closing(input, 5, tol).unwrap());
            assert!(
                closed == closing(input, 5, tol).unwrap(),
                "seed {seed} closing"
            );
        }
    }
    let s = memo.stats();
    assert!(s.hits > s.misses / 4, "{s:?}");
}

#[test]
fn memo_keeps_to_its_size() {
    // Too small for most of what it is given: nothing wrong may come of what is dropped.
    let memo = Memo::new(6_000);
    let tol = ArcTol::new(2, Side::Outside);
    for seed in 0..4u64 {
        let rings = scene(seed, 300, 130, 9);
        let got = memo.opening(&rings, 4, tol).unwrap();
        assert!(got == opening(&rings, 4, tol).unwrap(), "seed {seed}");
        assert!(memo.stats().bytes <= 6_000, "{:?}", memo.stats());
    }
    assert!(memo.stats().entries > 0);
    memo.clear();
    assert_eq!(memo.stats(), MemoStats::default());
    // A memo that can hold nothing is still a memo.
    let none = Memo::new(0);
    let rings = scene(9, 300, 130, 9);
    assert!(none.opening(&rings, 4, tol).unwrap() == opening(&rings, 4, tol).unwrap());
    assert_eq!(none.stats().entries, 0);
}

#[test]
fn memo_shared_between_threads() {
    fn assert_send_sync<T: Send + Sync + Clone>() {}
    assert_send_sync::<Memo>();
    let memo = Memo::default();
    let tol = ArcTol::new(2, Side::Inside);
    let scenes: Vec<Vec<Ring>> = (0..4).map(|k| scene(k % 2, 300, 130, 9)).collect();
    let want: Vec<PolygonSet> = scenes.iter().map(|s| opening(s, 5, tol).unwrap()).collect();
    std::thread::scope(|sc| {
        for (s, w) in scenes.iter().zip(&want) {
            let memo = memo.clone();
            sc.spawn(move || {
                for _ in 0..3 {
                    assert!(memo.opening(s, 5, tol).unwrap() == *w);
                }
            });
        }
    });
    assert!(memo.stats().hits > 0);
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(48))]

    /// A cluster's result moves with it: rings moved as one give the result moved.
    #[test]
    fn booleans_move_with_their_input(
        seed in 0u64..1000,
        dx in -(1i64 << 38)..(1i64 << 38),
        dy in -(1i64 << 38)..(1i64 << 38),
        k in 0u8..16,
    ) {
        let a = scene(seed, 40, 61, 5);
        let b = scene(seed + 7, 40, 61, 4);
        let shift = |rs: &[Ring]| -> Vec<Ring> {
            rs.iter()
                .map(|r| r.0.iter().map(|p| Point::new(p.x + dx, p.y + dy)).collect())
                .collect()
        };
        let run = |a: &[Ring], b: &[Ring]| {
            Boolean::new()
                .subject(a, fill(k))
                .clip(b, fill(k / 4))
                .op(op(k))
                .execute()
                .unwrap()
        };
        let here = run(&a, &b);
        let there = run(&shift(&a), &shift(&b));
        let moved: PolygonSet = here
            .iter()
            .map(|p| {
                let ring = |r: &Ring| -> Ring {
                    r.0.iter().map(|q| Point::new(q.x + dx, q.y + dy)).collect()
                };
                Polygon::new(ring(&p.outer), p.holes.iter().map(ring).collect())
            })
            .collect();
        prop_assert!(there == moved);
    }

    /// Scenes of repeated shapes through a memo of any size, twice, against one piece.
    #[test]
    fn memo_matches_one_piece(
        seed in 0u64..1000,
        n in 150usize..400,
        size in prop_oneof![Just(0usize), Just(4_000), Just(1 << 22)],
        k in 0u8..16,
    ) {
        let memo = Memo::new(size);
        let subject = scene(seed, n, 97, 7);
        let clip = scene(seed + 1, n / 2, 113, 6);
        let run = |memo: Option<&Memo>| {
            let mut b = Boolean::new()
                .subject(&subject, fill(k))
                .clip(&clip, fill(k / 4))
                .op(op(k))
                .monolithic(memo.is_none());
            b.set_memo(memo);
            b.execute_tree().unwrap()
        };
        let plain = run(None);
        prop_assert!(verified(|| run(Some(&memo))) == plain);
        prop_assert!(verified(|| run(Some(&memo))) == plain);
    }
}
