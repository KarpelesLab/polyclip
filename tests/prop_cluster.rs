//! Booleans computed cluster by cluster (independent groups of rings) must be identical to
//! the one-piece computation, for every op, fill rule and input.

use polyclip::*;
use proptest::prelude::*;

#[path = "common/corpus.rs"]
mod corpus;

/// One ring of a scene, placed in a grid cell so that neighbours are disjoint, near (within
/// a few units), touching, overlapping or nested.
#[derive(Clone, Debug)]
enum Shape {
    Rect {
        w: i64,
        h: i64,
    },
    /// Random vertices within the cell (often self-intersecting).
    Poly(Vec<(i64, i64)>),
    /// A rectangle around a block of cells (encloses other rings).
    Frame {
        cells: i64,
        margin: i64,
    },
    /// A thin sliver, diagonal or not.
    Sliver {
        dx: i64,
        dy: i64,
    },
    /// Concentric squares `gap` apart (separate clusters unless the gap is tiny), each
    /// reversed or not.
    Nest {
        depth: usize,
        gap: i64,
        rev: u32,
    },
}

#[derive(Clone, Debug)]
struct Item {
    cell: (i64, i64),
    jitter: (i64, i64),
    shape: Shape,
    reverse: bool,
    clip: bool,
    tag: u64,
}

fn item(pitch: i64) -> impl Strategy<Value = Item> {
    let shape = prop_oneof![
        4 => (1..pitch + 4, 1..pitch + 4).prop_map(|(w, h)| Shape::Rect { w, h }),
        3 => prop::collection::vec((0..pitch + 3, 0..pitch + 3), 3..7).prop_map(Shape::Poly),
        1 => (1i64..4, -3i64..4).prop_map(|(cells, margin)| Shape::Frame { cells, margin }),
        1 => (-pitch..pitch, -pitch..pitch).prop_map(|(dx, dy)| Shape::Sliver { dx, dy }),
        1 => (2usize..6, 0i64..4, any::<u32>())
            .prop_map(|(depth, gap, rev)| Shape::Nest { depth, gap, rev }),
    ];
    (
        (0i64..5, 0i64..5),
        (-3i64..4, -3i64..4),
        shape,
        any::<bool>(),
        any::<bool>(),
        0u64..4,
    )
        .prop_map(|(cell, jitter, shape, reverse, clip, tag)| Item {
            cell,
            jitter,
            shape,
            reverse,
            clip,
            tag,
        })
}

fn rings_of(it: &Item, pitch: i64) -> Vec<TaggedRing> {
    let Shape::Nest { depth, gap, rev } = it.shape else {
        return vec![ring_of(it, pitch)];
    };
    let (cx, cy) = (
        it.cell.0 * pitch + it.jitter.0,
        it.cell.1 * pitch + it.jitter.1,
    );
    (0..depth)
        .map(|k| {
            let r = 2 + k as i64 * (gap + 1);
            let mut pts = vec![
                Point::new(cx - r, cy - r),
                Point::new(cx + r, cy - r),
                Point::new(cx + r, cy + r),
                Point::new(cx - r, cy + r),
            ];
            if rev >> k & 1 != 0 {
                pts.reverse();
            }
            TaggedRing {
                tags: vec![it.tag * 16 + k as u64; 4],
                points: pts,
            }
        })
        .collect()
}

fn ring_of(it: &Item, pitch: i64) -> TaggedRing {
    let (x0, y0) = (
        it.cell.0 * pitch + it.jitter.0,
        it.cell.1 * pitch + it.jitter.1,
    );
    let mut pts: Vec<Point> = match &it.shape {
        Shape::Rect { w, h } => vec![
            Point::new(x0, y0),
            Point::new(x0 + w, y0),
            Point::new(x0 + w, y0 + h),
            Point::new(x0, y0 + h),
        ],
        Shape::Poly(v) => v.iter().map(|&(x, y)| Point::new(x0 + x, y0 + y)).collect(),
        Shape::Frame { cells, margin } => {
            let (a, b) = (x0 - margin, y0 - margin);
            let s = cells * pitch + 2 * margin;
            vec![
                Point::new(a, b),
                Point::new(a + s, b),
                Point::new(a + s, b + s),
                Point::new(a, b + s),
            ]
        }
        Shape::Sliver { dx, dy } => vec![
            Point::new(x0, y0),
            Point::new(x0 + dx, y0 + dy),
            Point::new(x0 + dx + 1, y0 + dy),
            Point::new(x0 + 1, y0 + 1),
        ],
        Shape::Nest { .. } => unreachable!(),
    };
    if it.reverse {
        pts.reverse();
    }
    // Tags vary along the ring so that provenance is checked too.
    let n = pts.len() as u64;
    let tags = (0..n).map(|k| it.tag * 16 + k % 3).collect();
    TaggedRing { points: pts, tags }
}

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

/// Runs the operation clustered (forced) and in one piece; both trees must be identical.
fn same(
    items: &[Item],
    pitch: i64,
    scale: i64,
    o: Op,
    ra: FillRule,
    rb: FillRule,
    keep: bool,
) -> core::result::Result<(), TestCaseError> {
    let mut subj: Vec<TaggedRing> = Vec::new();
    let mut clip: Vec<TaggedRing> = Vec::new();
    for it in items {
        for mut r in rings_of(it, pitch) {
            for p in r.points.iter_mut() {
                *p = Point::new(p.x * scale, p.y * scale);
            }
            if it.clip { clip.push(r) } else { subj.push(r) }
        }
    }
    let run = |force: bool, mono: bool| {
        Boolean::new()
            .subject(&subj, ra)
            .clip(&clip, rb)
            .op(o)
            .keep_collinear(keep)
            .monolithic(mono)
            .force_clusters(force)
            .execute_tree()
            .unwrap()
    };
    let one = run(false, true);
    let clustered = run(true, false);
    prop_assert_eq!(&clustered, &one);
    // The automatic choice too.
    prop_assert_eq!(&run(false, false), &one);
    Ok(())
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(2000))]

    /// Small pitch: rings often touch, overlap or come within a unit or two.
    #[test]
    fn clustered_dense(items in prop::collection::vec(item(6), 1..14), o in op(),
                       ra in rule(), rb in rule(), keep in any::<bool>()) {
        same(&items, 6, 1, o, ra, rb, keep)?;
    }

    /// Larger pitch: mostly separate rings, some nested in frames.
    #[test]
    fn clustered_sparse(items in prop::collection::vec(item(12), 1..20), o in op(),
                        ra in rule(), rb in rule(), keep in any::<bool>()) {
        same(&items, 12, 1, o, ra, rb, keep)?;
    }

    /// Scaled up: no snap rounding interactions between cells, many lone rings.
    #[test]
    fn clustered_scaled(items in prop::collection::vec(item(10), 1..20), o in op(),
                        ra in rule(), rb in rule(), keep in any::<bool>()) {
        same(&items, 10, 1000, o, ra, rb, keep)?;
    }
}

/// Lone rings at exactly the clustering margin: gaps of 0 to 8 units, along and across the
/// axes and diagonally, with crossings whose rounding lands near the neighbour.
#[test]
fn near_margin() {
    let p = Point::new;
    for gap in 0..9 {
        for o in [Op::Union, Op::Intersection, Op::Difference, Op::Xor] {
            let a = Ring::from([(0, 0), (10, 0), (10, 10), (0, 10)]);
            // A self-crossing bow tie near `a` (its crossing is rounded).
            let bow: Ring = [
                p(10 + gap, 0),
                p(17 + gap, 7),
                p(17 + gap, 0),
                p(10 + gap, 7),
            ]
            .into_iter()
            .collect();
            let diag: Ring = [
                p(10 + gap, 10 + gap),
                p(20 + gap, 13 + gap),
                p(13 + gap, 20 + gap),
            ]
            .into_iter()
            .collect();
            let thin: Ring = [p(-gap - 1, 0), p(-gap, 0), p(-gap, 30), p(-gap - 1, 31)]
                .into_iter()
                .collect();
            let subj = vec![a.clone(), bow.clone(), thin];
            let clip = vec![diag, a];
            let run = |force: bool| {
                Boolean::new()
                    .subject(&subj, FillRule::NonZero)
                    .clip(&clip, FillRule::EvenOdd)
                    .op(o)
                    .monolithic(!force)
                    .force_clusters(force)
                    .execute_tree()
                    .unwrap()
            };
            assert_eq!(run(true), run(false), "gap {gap}, {o:?}");
        }
    }
}

/// The real corpus: one polygon with 2153 holes, with small rectangles added (some crossing
/// holes, some inside holes, some in copper), for every op.
#[test]
fn corpus_with_spokes() {
    let fill = corpus::load("cadlab_gnd_in1.pclp");
    let b = fill[0].outer.bbox().unwrap();
    let mut s = 7u64;
    let mut rnd = |m: i64| {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((s >> 33) as i64).rem_euclid(m)
    };
    let spokes: Vec<Ring> = (0..60)
        .map(|_| {
            let x = b.min.x + rnd(b.width());
            let y = b.min.y + rnd(b.height());
            let (w, h) = (1 + rnd(2_000_000), 1 + rnd(2_000_000));
            Ring::from([(x, y), (x + w, y), (x + w, y + h), (x, y + h)])
        })
        .collect();
    for o in [Op::Union, Op::Intersection, Op::Difference, Op::Xor] {
        let run = |mono: bool| {
            Boolean::new()
                .subject(&fill, FillRule::NonZero)
                .clip(&spokes, FillRule::NonZero)
                .op(o)
                .monolithic(mono)
                .execute_tree()
                .unwrap()
        };
        assert_eq!(run(false), run(true), "{o:?}");
    }
    let norm = |mono: bool| {
        Boolean::new()
            .subject(&fill, FillRule::NonZero)
            .monolithic(mono)
            .execute_tree()
            .unwrap()
    };
    assert_eq!(norm(false), norm(true));
}

/// Larger random boards (above the size where clustering kicks in on its own): a zone
/// minus or plus many small polygons, some overlapping, some touching, some nested.
#[test]
fn random_boards() {
    let mut s = 0xb0a2du64;
    let mut rnd = move |m: i64| {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((s >> 33) as i64).rem_euclid(m)
    };
    for it in 0..8 {
        let field = 2_000 + rnd(40_000);
        let zone = Ring::from([(0, 0), (field, 0), (field, field), (0, field)]);
        let mut obst: Vec<TaggedRing> = Vec::new();
        for k in 0..(150 + rnd(250)) {
            let (cx, cy) = (rnd(field), rnd(field));
            let r = (3 + rnd(300)) as f64;
            let n = 3 + rnd(40) as usize;
            let mut pts: Vec<Point> = (0..n)
                .map(|j| {
                    let a = std::f64::consts::TAU * j as f64 / n as f64;
                    Point::new(
                        cx + (r * a.cos()).round() as i64,
                        cy + (r * a.sin()).round() as i64,
                    )
                })
                .collect();
            if rnd(4) == 0 {
                pts.reverse();
            }
            if rnd(8) == 0 {
                // A bow tie or spike now and then.
                pts.swap(0, n / 2);
            }
            obst.push(TaggedRing::uniform(Ring(pts), k as u64 % 5));
        }
        let rules = [
            FillRule::NonZero,
            FillRule::EvenOdd,
            FillRule::Positive,
            FillRule::Negative,
        ];
        for o in [Op::Union, Op::Intersection, Op::Difference, Op::Xor] {
            let rb = rules[rnd(4) as usize];
            let run = |mono: bool| {
                Boolean::new()
                    .subject(&zone, FillRule::NonZero)
                    .clip(&obst, rb)
                    .op(o)
                    .monolithic(mono)
                    .execute_tree()
                    .unwrap()
            };
            assert_eq!(run(false), run(true), "board {it}, {o:?}, {rb:?}");
        }
    }
}
