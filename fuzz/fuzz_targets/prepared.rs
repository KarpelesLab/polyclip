//! `Prepared` geometries: every query (`locate`, `intersects`, `contains`,
//! `distance_less_than`, `distance` with its closest points) must give exactly what the free
//! function gives, on raw (possibly invalid) and canonical geometry of every kind. The
//! prepared geometry is also tiled (copies side by side) so that the indexed paths, not
//! only the small-input fallbacks, are taken.
#![no_main]

use libfuzzer_sys::fuzz_target;
use polyclip::*;
use polyclip_fuzz::Gen;

fn same_closest(a: Option<Closest>, b: Option<Closest>) -> bool {
    match (a, b) {
        (None, None) => true,
        (Some(a), Some(b)) => {
            a.sq == b.sq
                && a.a.x.to_bits() == b.a.x.to_bits()
                && a.a.y.to_bits() == b.a.y.to_bits()
                && a.b.x.to_bits() == b.b.x.to_bits()
                && a.b.y.to_bits() == b.b.y.to_bits()
        }
        _ => false,
    }
}

fn check<A: Preparable + ?Sized, B: Geometry + ?Sized>(
    p: &Prepared<'_, A>,
    b: &B,
    pts: &[Point],
    d: i64,
) {
    let a = p.geometry();
    for &q in pts {
        assert_eq!(p.locate(q), locate(a, q), "locate {q:?}");
    }
    assert_eq!(p.intersects(b), intersects(a, b), "intersects");
    assert_eq!(p.contains(b), contains(a, b), "contains");
    assert_eq!(
        p.distance_less_than(b, d),
        distance_less_than(a, b, d),
        "distance_less_than {d}"
    );
    assert!(same_closest(p.distance(b), distance(a, b)), "distance");
}

fn all<A: Preparable + ?Sized>(a: &A, bs: &Bs, pts: &[Point], d: i64) {
    let p = Prepared::new(a);
    check(&p, &bs.ring, pts, d);
    check(&p, &bs.polys, pts, d);
    check(&p, &bs.path, pts, d);
    check(&p, &bs.seg, pts, d);
    check(&p, &bs.point, pts, d);
}

struct Bs {
    ring: Ring,
    polys: Vec<Polygon>,
    path: Path,
    seg: Segment,
    point: Point,
}

/// `k` copies of the rings side by side, `step` apart along x (when that stays in range).
fn tile(rings: &[Ring], k: i64, step: i64) -> Vec<Ring> {
    let ok = |p: &Point| p.x.unsigned_abs() as i128 + (k * step) as i128 <= MAX_COORD as i128;
    if !rings.iter().all(|r| r.iter().all(ok)) {
        return rings.to_vec();
    }
    (0..k)
        .flat_map(|i| {
            rings
                .iter()
                .map(move |r| r.iter().map(|p| Point::new(p.x + i * step, p.y)).collect())
        })
        .collect()
}

fuzz_target!(|data: &[u8]| {
    let mut g = Gen::new(data);
    let ra = g.rings(4, 12);
    let rb = g.rings(2, 6);
    let path = g.path(6);
    let p = g.point();
    let q = g.point();
    let polys = g.raw_polygons(2);
    let d = g.delta().unsigned_abs().min(MAX_COORD as u64) as i64;
    polyclip_fuzz::dump(
        "(ra, rb, path, p, q, polys, d)",
        &(&ra, &rb, &path, p, q, &polys, d),
    );
    let bs = Bs {
        ring: rb.first().cloned().unwrap_or_default(),
        polys: rb.iter().cloned().map(Polygon::from).collect(),
        path: path.clone(),
        seg: Segment::new(p, q),
        point: q,
    };
    let mut pts = vec![p, q];
    pts.extend(path.iter().copied());
    pts.extend(ra.iter().flat_map(|r| r.iter().copied()).take(16));

    // Raw geometry of every kind.
    for r in &ra {
        all(r, &bs, &pts, d);
    }
    all(&polys, &bs, &pts, d);
    for poly in &polys {
        all(poly, &bs, &pts, d);
    }
    all(&path, &bs, &pts, d);
    all(&Segment::new(p, q), &bs, &pts, d);
    all(&p, &bs, &pts, d);
    let tree = PolyTree {
        nodes: ra
            .iter()
            .enumerate()
            .map(|(i, r)| PolyNode {
                ring: r.clone(),
                tags: vec![0; r.len()],
                is_hole: i % 2 == 1,
                parent: None,
                children: Vec::new(),
            })
            .collect(),
        roots: Vec::new(),
    };
    all(&tree, &bs, &pts, d);

    // Tiled: enough segments for the indexed searches.
    let n: usize = ra.iter().map(|r| r.len()).sum();
    if n > 0 {
        let k = (600 / n + 1).min(200) as i64;
        let step = 2 * g.extent() + 3;
        let big = tile(&ra, k, step);
        let set: Vec<Polygon> = big.iter().cloned().map(Polygon::from).collect();
        all(&set, &bs, &pts, d);
        let poly = Polygon::new(big[0].clone(), big[1..].to_vec());
        all(&poly, &bs, &pts, d);
        if let Ok(c) = union_all(&big, g.fill_rule()) {
            all(&c, &bs, &pts, d);
        }
    }
    if let Ok(c) = union_all(&ra, g.fill_rule()) {
        all(&c, &bs, &pts, d);
    }
});
