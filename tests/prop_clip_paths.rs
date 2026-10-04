//! Open-path clipping computed cluster by cluster (only the clip rings near the paths) must
//! be identical to the one-piece computation.

use polyclip::*;
use proptest::prelude::*;

#[path = "common/corpus.rs"]
mod corpus;

const PITCH: i64 = 10;

/// A clip ring: a few vertices around a grid cell (neighbours near, touching, overlapping,
/// self-intersecting), or a rectangle around a block of cells (enclosing others).
fn ring() -> impl Strategy<Value = Ring> {
    prop_oneof![
        8 => ((0i64..40, 0i64..40), prop::collection::vec((-2i64..PITCH + 3, -2i64..PITCH + 3), 4..9))
            .prop_map(|((cx, cy), v)| v.into_iter().map(|(x, y)| Point::new(cx * PITCH + x, cy * PITCH + y)).collect()),
        3 => ((0i64..40, 0i64..40), 1..PITCH, 1..PITCH, any::<bool>()).prop_map(|((cx, cy), w, h, rev)| {
            let (x, y) = (cx * PITCH, cy * PITCH);
            let mut r = Ring::from([(x, y), (x + w, y), (x + w, y + h), (x, y + h)]);
            if rev { r.0.reverse() }
            r
        }),
        1 => ((0i64..40, 0i64..40), 1i64..8, -2i64..3).prop_map(|((cx, cy), n, m)| {
            let (x, y, s) = (cx * PITCH - m, cy * PITCH - m, n * PITCH + 2 * m);
            Ring::from([(x, y), (x + s, y), (x + s, y + s), (x, y + s)])
        }),
    ]
}

fn path() -> impl Strategy<Value = TaggedPath> {
    let short = (
        (0i64..400, 0i64..400),
        prop::collection::vec((-6i64..7, -6i64..7), 1..5),
    )
        .prop_map(|((x, y), steps)| {
            let mut p = Point::new(x, y);
            let mut pts = vec![p];
            for (dx, dy) in steps {
                p = Point::new(p.x + dx, p.y + dy);
                pts.push(p);
            }
            pts
        });
    let long = prop::collection::vec((-5i64..405, -5i64..405), 1..5)
        .prop_map(|v| v.into_iter().map(Point::from).collect::<Vec<_>>());
    (prop_oneof![3 => short, 1 => long], 0u64..5).prop_map(|(points, t)| TaggedPath {
        tags: (0..points.len().saturating_sub(1) as u64)
            .map(|k| t * 8 + k)
            .collect(),
        points,
    })
}

fn rule() -> impl Strategy<Value = FillRule> {
    prop_oneof![
        Just(FillRule::EvenOdd),
        Just(FillRule::NonZero),
        Just(FillRule::Positive),
        Just(FillRule::Negative)
    ]
}

fn one_piece<P: PathSource + ?Sized, C: RingSource + ?Sized>(
    p: &P,
    c: &C,
    r: FillRule,
) -> ClippedPaths {
    set_always_monolithic(true);
    let out = clip_paths(p, c, r);
    set_always_monolithic(false);
    out.unwrap()
}

/// 48 cases unless `PROPTEST_CASES` says otherwise.
fn cases() -> u32 {
    std::env::var("PROPTEST_CASES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(48)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(cases()))]

    #[test]
    fn clustered_clip_paths(rings in prop::collection::vec(ring(), 200..320), paths in prop::collection::vec(path(), 1..24), r in rule(), scale in prop_oneof![Just(1i64), Just(1000)]) {
        let rings: Vec<Ring> = rings
            .into_iter()
            .map(|g| g.0.into_iter().map(|p| Point::new(p.x * scale, p.y * scale)).collect())
            .collect();
        let paths: Vec<TaggedPath> = paths
            .into_iter()
            .map(|mut p| {
                for q in p.points.iter_mut() {
                    *q = Point::new(q.x * scale, q.y * scale);
                }
                p
            })
            .collect();
        let want = one_piece(&paths, &rings, r);
        prop_assert_eq!(&clip_paths(&paths, &rings, r).unwrap(), &want);
        // Each path alone, and the clip as canonical polygons.
        for p in paths.iter().take(4) {
            prop_assert_eq!(clip_paths(p, &rings, r).unwrap(), one_piece(p, &rings, r));
        }
        let polys = union_all(&rings, r).unwrap();
        prop_assert_eq!(clip_paths(&paths, &polys, r).unwrap(), one_piece(&paths, &polys, r));
    }
}

#[test]
fn clustered_clip_paths_corpus() {
    let fill = corpus::load("cadlab_gnd_in1.pclp");
    let b = fill.bbox().unwrap();
    let mut s = 3u64;
    let mut rnd = |m: i64| {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((s >> 33) as i64).rem_euclid(m)
    };
    let tracks: Vec<Path> = (0..60)
        .map(|_| {
            let mut p = Point::new(b.min.x + rnd(b.width()), b.min.y + rnd(b.height()));
            let mut pts = vec![p];
            for _ in 0..1 + rnd(4) {
                let len = 200_000 + rnd(4_000_000);
                let (dx, dy) = [
                    (1, 0),
                    (1, 1),
                    (0, 1),
                    (-1, 1),
                    (-1, 0),
                    (-1, -1),
                    (0, -1),
                    (1, -1),
                ][rnd(8) as usize];
                p = Point::new(p.x + dx * len, p.y + dy * len);
                pts.push(p);
            }
            Path(pts)
        })
        .collect();
    assert_eq!(
        clip_paths(&tracks, &fill, FillRule::NonZero).unwrap(),
        one_piece(&tracks, &fill, FillRule::NonZero)
    );
    for t in tracks.iter().take(4) {
        assert_eq!(
            clip_paths(t, &fill, FillRule::EvenOdd).unwrap(),
            one_piece(t, &fill, FillRule::EvenOdd)
        );
    }
}
