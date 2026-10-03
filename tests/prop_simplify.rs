use polyclip::*;
use proptest::prelude::*;
use std::collections::HashMap;

/// Arbitrary (usually self-intersecting) rings.
fn rings(range: i64, max_rings: usize, max_pts: usize) -> impl Strategy<Value = Vec<Ring>> {
    prop::collection::vec(
        prop::collection::vec((-range..=range, -range..=range), 3..max_pts)
            .prop_map(|v| v.into_iter().map(Point::from).collect::<Ring>()),
        1..max_rings,
    )
}

/// A star-shaped ring around `(cx, cy)` with `n` vertices and per-vertex radial noise.
fn wobbly(cx: i64, cy: i64, r: f64, noise: &[f64]) -> Ring {
    let n = noise.len();
    (0..n)
        .map(|k| {
            let a = k as f64 * std::f64::consts::TAU / n as f64;
            let rr = r * (1.0 + noise[k]);
            Point::new(
                cx + (rr * a.cos()).round() as i64,
                cy + (rr * a.sin()).round() as i64,
            )
        })
        .collect()
}

fn wobbly_circles(max: usize) -> impl Strategy<Value = Vec<Ring>> {
    prop::collection::vec(
        (
            0i64..3000,
            0i64..3000,
            20.0f64..800.0,
            0.0f64..0.3,
            prop::collection::vec(-1.0f64..1.0, 8..300),
        )
            .prop_map(|(cx, cy, r, amp, noise)| {
                let noise: Vec<f64> = noise.iter().map(|v| v * amp).collect();
                wobbly(cx, cy, r, &noise)
            }),
        1..max,
    )
}

fn perimeter(set: &[Polygon]) -> f64 {
    set.iter()
        .flat_map(|p| p.rings())
        .flat_map(|r| r.edges())
        .map(|(a, b)| ((b.x - a.x) as f64).hypot((b.y - a.y) as f64))
        .sum()
}

fn seg_dist(p: Point, a: Point, b: Point) -> f64 {
    let (px, py) = ((p.x - a.x) as f64, (p.y - a.y) as f64);
    let (dx, dy) = ((b.x - a.x) as f64, (b.y - a.y) as f64);
    let l2 = dx * dx + dy * dy;
    let t = if l2 == 0.0 {
        0.0
    } else {
        ((px * dx + py * dy) / l2).clamp(0.0, 1.0)
    };
    (px - t * dx).hypot(py - t * dy)
}

/// Checks that `out` is `ring` with some vertices removed (same cyclic order) and every
/// removed vertex within `tol` of its replacing edge.
fn is_simplification(ring: &Ring, out: &Ring, tol: i64) -> std::result::Result<(), String> {
    let n = ring.len();
    let start = ring
        .iter()
        .position(|&p| p == out[0])
        .ok_or("start not found")?;
    let mut k = 0;
    let mut matched = Vec::new();
    for t in 0..n {
        let i = (start + t) % n;
        if k < out.len() && ring[i] == out[k] {
            matched.push(i);
            k += 1;
        }
    }
    if k != out.len() {
        return Err("not a subsequence".into());
    }
    let slack = 1e-6 * (1.0 + tol as f64);
    for w in 0..matched.len() {
        let (i, j) = (matched[w], matched[(w + 1) % matched.len()]);
        let (a, b) = (ring[i], ring[j]);
        let mut x = (i + 1) % n;
        while x != j {
            let d = seg_dist(ring[x], a, b);
            if d > tol as f64 + slack {
                return Err(format!("vertex {:?} at {} > {}", ring[x], d, tol));
            }
            x = (x + 1) % n;
        }
    }
    Ok(())
}

fn check(input: &[Polygon], tol: i64) -> std::result::Result<PolygonSet, TestCaseError> {
    let out = simplify_polygons(input, tol);
    prop_assert_eq!(check_canonical(&out, true), Ok(()));
    prop_assert_eq!(out.len(), input.len());
    let mut hole_counts: Vec<usize> = input.iter().map(|p| p.holes.len()).collect();
    let mut out_counts: Vec<usize> = out.iter().map(|p| p.holes.len()).collect();
    hole_counts.sort_unstable();
    out_counts.sort_unstable();
    prop_assert_eq!(hole_counts, out_counts);
    // Every output ring is a simplification of some input ring with the same role.
    let mut by_point: HashMap<(Point, bool), Vec<&Ring>> = HashMap::new();
    for poly in input {
        for (k, r) in poly.rings().enumerate() {
            for &p in r.iter() {
                by_point.entry((p, k == 0)).or_default().push(r);
            }
        }
    }
    for poly in &out {
        for (k, r) in poly.rings().enumerate() {
            prop_assert!(r.len() >= 3);
            let cands = by_point.get(&(r[0], k == 0)).cloned().unwrap_or_default();
            let errs: Vec<String> = cands
                .iter()
                .filter_map(|c| is_simplification(c, r, tol).err())
                .collect();
            prop_assert!(errs.len() < cands.len(), "ring {:?}: {:?}", r, errs);
        }
    }
    // Area change bounded by the swept stadiums.
    let a_in: i128 = input.iter().map(|p| p.signed_area2()).sum();
    let a_out: i128 = out.iter().map(|p| p.signed_area2()).sum();
    let nv: usize = input.iter().map(|p| p.vertex_count()).sum();
    let t = tol.max(0) as f64;
    let bound = 2.0 * (2.0 * t * perimeter(input) + 4.0 * t * t * nv as f64) + 1.0;
    prop_assert!(
        ((a_in - a_out) as f64).abs() <= bound,
        "area {} -> {} bound {}",
        a_in,
        a_out,
        bound
    );
    // Tolerance zero only drops collinear vertices: canonical input is unchanged.
    if tol <= 0 {
        prop_assert_eq!(&out, &input.to_vec());
    }
    Ok(out)
}

proptest! {
    #[test]
    fn random_rings(rs in rings(200, 6, 9), tol in -2i64..80) {
        let set = union_all(&rs, FillRule::NonZero).unwrap();
        let out = check(&set, tol)?;
        // Idempotent enough: simplifying again stays valid.
        check(&out, tol)?;
    }

    #[test]
    fn random_rings_small(rs in rings(15, 8, 9), tol in 0i64..30) {
        // Tiny coordinate range: many pinch points, touching rings and holes.
        let set = union_all(&rs, FillRule::EvenOdd).unwrap();
        check(&set, tol)?;
    }

    #[test]
    fn random_rings_large(rs in rings(1 << 40, 5, 9), tol in prop_oneof![0i64..1000, 0i64..(1 << 42), Just(i64::MAX)]) {
        let set = union_all(&rs, FillRule::NonZero).unwrap();
        check(&set, tol)?;
    }

    #[test]
    fn wobbly_sets(a in wobbly_circles(6), b in wobbly_circles(5), tol in prop_oneof![0i64..10, 0i64..200, 0i64..3000]) {
        // Circles minus circles: holes, islands and pinch points.
        let sa = union_all(&a, FillRule::NonZero).unwrap();
        let sb = union_all(&b, FillRule::NonZero).unwrap();
        let set = boolean(Op::Difference, &sa, &sb, FillRule::NonZero).unwrap();
        check(&set, tol)?;
        let x = boolean(Op::Xor, &sa, &sb, FillRule::NonZero).unwrap();
        check(&x, tol)?;
    }

    #[test]
    fn path_props(pts in prop::collection::vec((-1000i64..1000, -1000i64..1000), 0..40), tol in -2i64..300) {
        let path: Path = pts.into_iter().map(Point::from).collect();
        let out = simplify_path(&path, tol);
        let mut dedup = path.to_vec();
        dedup.dedup();
        if !dedup.is_empty() {
            prop_assert_eq!(out.first(), dedup.first());
            prop_assert_eq!(out.last(), dedup.last());
        }
        // Subsequence with the Hausdorff bound.
        let mut k = 0;
        let mut last = 0;
        let tolf = tol.max(0) as f64 + 1e-6;
        for (i, &p) in dedup.iter().enumerate() {
            if k < out.len() && p == out[k] {
                if k > 0 {
                    for &q in &dedup[last + 1..i] {
                        prop_assert!(seg_dist(q, out[k - 1], p) <= tolf);
                    }
                }
                last = i;
                k += 1;
            }
        }
        prop_assert_eq!(k, out.len());
    }

    #[test]
    fn garbage_never_panics(rs in rings(50, 4, 12), holes in rings(50, 3, 6), tol in any::<i64>()) {
        let polys: PolygonSet = rs.into_iter().map(|r| Polygon::new(r, holes.clone())).collect();
        let out = simplify_polygons(&polys, tol);
        prop_assert_eq!(out.len(), polys.len());
    }
}

/// 100k+ vertices, 200 rings: simplification must be fast (timed in release builds only).
#[test]
fn perf_large_set() {
    let mut rs = Vec::new();
    let mut seed = 12345u64;
    let mut rnd = || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((seed >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
    };
    for gx in 0..10 {
        for gy in 0..10 {
            let (cx, cy) = (gx * 100_000, gy * 100_000);
            let outer: Vec<f64> = (0..700).map(|_| 0.02 * rnd()).collect();
            rs.push(wobbly(cx, cy, 45_000.0, &outer));
            let inner: Vec<f64> = (0..350).map(|_| 0.02 * rnd()).collect();
            let mut h = wobbly(cx, cy, 20_000.0, &inner);
            h.reverse_orientation();
            rs.push(h);
        }
    }
    let set = union_all(&rs, FillRule::NonZero).unwrap();
    let nv: usize = set.iter().map(|p| p.vertex_count()).sum();
    assert!(nv > 100_000, "{nv}");
    for tol in [0, 100, 1000, 100_000] {
        let t = std::time::Instant::now();
        let out = simplify_polygons(&set, tol);
        let el = t.elapsed();
        let ov: usize = out.iter().map(|p| p.vertex_count()).sum();
        eprintln!("tol {tol}: {nv} -> {ov} vertices in {el:?}");
        if !cfg!(debug_assertions) {
            assert!(el.as_secs_f64() < 2.0, "too slow: {el:?}");
        }
        assert_eq!(validate_set(&out), Ok(()));
        assert_eq!(out.len(), set.len());
    }
}
