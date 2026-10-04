//! The parallel boolean engine (banded sweep, dense snap rounding) must give exactly the
//! same result as the original sequential implementation, for every op, fill rule and
//! input, however finely the work is split.

use polyclip::*;
use proptest::prelude::*;

#[path = "common/corpus.rs"]
mod corpus;

fn rule() -> impl Strategy<Value = FillRule> {
    prop_oneof![
        Just(FillRule::EvenOdd),
        Just(FillRule::NonZero),
        Just(FillRule::Positive),
        Just(FillRule::Negative)
    ]
}

fn op() -> impl Strategy<Value = Op> {
    prop_oneof![
        Just(Op::Union),
        Just(Op::Intersection),
        Just(Op::Difference),
        Just(Op::Xor)
    ]
}

/// A ring with vertices in `[0, range)` (often self-intersecting), or an axis-parallel
/// rectangle (coincident edges), with varying tags.
fn ring(range: i64) -> impl Strategy<Value = TaggedRing> {
    let poly = prop::collection::vec((0..range, 0..range), 3..9);
    let rect = (0..range, 0..range, 1..range, 1..range)
        .prop_map(|(x, y, w, h)| vec![(x, y), (x + w, y), (x + w, y + h), (x, y + h)]);
    (prop_oneof![3 => poly, 1 => rect], 0u64..3, any::<bool>()).prop_map(|(v, tag, rev)| {
        let mut points: Vec<Point> = v.into_iter().map(|(x, y)| Point::new(x, y)).collect();
        if rev {
            points.reverse();
        }
        let tags = (0..points.len() as u64).map(|k| tag * 4 + k % 2).collect();
        TaggedRing { points, tags }
    })
}

fn run(
    engine: Engine,
    subj: &[TaggedRing],
    clip: &[TaggedRing],
    o: Op,
    ra: FillRule,
    rb: FillRule,
    keep: bool,
) -> PolyTree {
    set_engine(engine);
    let t = Boolean::new()
        .subject(subj, ra)
        .clip(clip, rb)
        .op(o)
        .keep_collinear(keep)
        .monolithic(true)
        .execute_tree()
        .unwrap();
    set_engine(Engine::Auto);
    t
}

fn same(
    subj: &[TaggedRing],
    clip: &[TaggedRing],
    o: Op,
    ra: FillRule,
    rb: FillRule,
    keep: bool,
) -> core::result::Result<(), TestCaseError> {
    let reference = run(Engine::Reference, subj, clip, o, ra, rb, keep);
    prop_assert_eq!(&run(Engine::Split, subj, clip, o, ra, rb, keep), &reference);
    prop_assert_eq!(&run(Engine::Auto, subj, clip, o, ra, rb, keep), &reference);
    Ok(())
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(1500))]

    /// Tiny coordinate range: crossings everywhere, heavy snap rounding.
    #[test]
    fn engine_tiny(subj in prop::collection::vec(ring(8), 0..6),
                   clip in prop::collection::vec(ring(8), 0..6),
                   o in op(), ra in rule(), rb in rule(), keep in any::<bool>()) {
        same(&subj, &clip, o, ra, rb, keep)?;
    }

    /// Small range: rounding interacts with nearby edges.
    #[test]
    fn engine_small(subj in prop::collection::vec(ring(40), 0..8),
                    clip in prop::collection::vec(ring(40), 0..8),
                    o in op(), ra in rule(), rb in rule(), keep in any::<bool>()) {
        same(&subj, &clip, o, ra, rb, keep)?;
    }

    /// Larger range: mostly exact crossings far apart.
    #[test]
    fn engine_wide(subj in prop::collection::vec(ring(100_000), 0..8),
                   clip in prop::collection::vec(ring(100_000), 0..8),
                   o in op(), ra in rule(), rb in rule(), keep in any::<bool>()) {
        same(&subj, &clip, o, ra, rb, keep)?;
    }
}

fn circle(cx: i64, cy: i64, r: f64, n: usize) -> Ring {
    (0..n)
        .map(|k| {
            let a = core::f64::consts::TAU * k as f64 / n as f64;
            Point::new(
                cx + (r * a.cos()).round() as i64,
                cy + (r * a.sin()).round() as i64,
            )
        })
        .collect()
}

/// Heavily overlapping circles (large enough for the automatic engine to split the work),
/// with every op against a second set.
#[test]
fn overlapping_circles() {
    let mut s = 42u64;
    let mut rnd = move |m: u64| {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (s >> 33) % m
    };
    for (count, field, r) in [(1500, 4_000_000, 150_000.0), (1500, 30_000, 900.0)] {
        let mut make = |count: usize| -> Vec<TaggedRing> {
            (0..count)
                .map(|k| {
                    let c = circle(rnd(field) as i64, rnd(field) as i64, r, 24);
                    TaggedRing::uniform(c, k as u64 % 7)
                })
                .collect()
        };
        let a = make(count);
        let b = make(count / 3);
        for (o, rb) in [
            (Op::Union, FillRule::NonZero),
            (Op::Difference, FillRule::EvenOdd),
            (Op::Xor, FillRule::Positive),
            (Op::Intersection, FillRule::NonZero),
        ] {
            let reference = run(Engine::Reference, &a, &b, o, FillRule::NonZero, rb, false);
            for e in [Engine::Split, Engine::Auto] {
                let t = run(e, &a, &b, o, FillRule::NonZero, rb, false);
                assert!(t == reference, "{e:?} {o:?} {rb:?} field {field}");
            }
        }
    }
}

/// The real corpus, normalized and with rectangles added, in every engine.
#[test]
fn corpus_engines() {
    let set = corpus::load("cadlab_gnd_in1.pclp");
    let b = set[0].bbox().unwrap();
    let fill: Vec<TaggedRing> = set
        .iter()
        .flat_map(|p| p.rings().map(|r| TaggedRing::uniform(r.clone(), 1)))
        .collect();
    let mut s = 9u64;
    let mut rnd = |m: i64| {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((s >> 33) as i64).rem_euclid(m)
    };
    let spokes: Vec<TaggedRing> = (0..200)
        .map(|_| {
            let x = b.min.x + rnd(b.width());
            let y = b.min.y + rnd(b.height());
            let (w, h) = (1 + rnd(3_000_000), 1 + rnd(3_000_000));
            TaggedRing::uniform(
                Ring::from([(x, y), (x + w, y), (x + w, y + h), (x, y + h)]),
                2,
            )
        })
        .collect();
    for o in [Op::Union, Op::Difference, Op::Xor] {
        let reference = run(
            Engine::Reference,
            &fill,
            &spokes,
            o,
            FillRule::NonZero,
            FillRule::NonZero,
            false,
        );
        for e in [Engine::Split, Engine::Auto] {
            let t = run(
                e,
                &fill,
                &spokes,
                o,
                FillRule::NonZero,
                FillRule::NonZero,
                false,
            );
            assert!(t == reference, "{e:?} {o:?}");
        }
    }
}
