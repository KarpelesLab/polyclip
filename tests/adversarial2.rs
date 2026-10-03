//! Second adversarial review.
//!
//! The first group reproduces bugs found by the review (all fixed; each documents expected vs
//! actual and the suspected root cause). The second group holds randomized differential
//! checks of code that was investigated and found correct (distance / intersection sweeps
//! against brute force, trapezoid coverage, triangulation and offset invariants).

use polyclip::*;
use std::cmp::Ordering;

fn p(x: i64, y: i64) -> Point {
    Point::new(x, y)
}

fn rect(x0: i64, y0: i64, w: i64, h: i64) -> Ring {
    Ring::from([(x0, y0), (x0 + w, y0), (x0 + w, y0 + h), (x0, y0 + h)])
}

/// Deterministic LCG for the self-contained randomized checks.
struct Rng(u64);

impl Rng {
    fn r(&mut self, m: i64) -> i64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 33) as i64).rem_euclid(m)
    }
}

// =======================================================================================
// Confirmed bugs
// =======================================================================================

/// Expected: a polygon whose hole has coordinates beyond `MAX_COORD` is out of range:
/// `in_range` is `false` and the infallible queries return their neutral result (`0`,
/// `None`, `false`), as the crate docs promise.
/// Actual: `in_range` returns `true`, and `area2`, `centroid`, `intersects`, `distance`
/// and `contains` compute with overflowing arithmetic (panics with "attempt to subtract /
/// multiply / add with overflow" in debug builds, garbage in release builds).
///
/// Root cause: `Polygon::bbox` (src/geom.rs:339, used by `impl Geometry for Polygon` in
/// src/query.rs:258) is the outer ring's box only, and every range check
/// (`in_range` src/query.rs:400, `area2` :676, `centroid` :689, `intersects` :557,
/// `distance` src/distance.rs:206, `contains` src/query.rs:607) checks just the bbox.
/// The holes are visited by `visit_segments` all the same. (Fallible operations check
/// every vertex and correctly report `CoordinateOutOfRange`.)
#[test]
fn out_of_range_hole_is_detected() {
    let big = i64::MAX;
    let q = Polygon::new(
        rect(0, 0, 10, 10),
        vec![Ring::from([(-big, -big), (big, -big), (0, big)])],
    );
    assert!(validate(&q).is_err()); // fallible paths do see the hole
    assert!(!in_range(&q)); // actual: true
    assert_eq!(area2(&q), 0); // actual: overflow panic (debug)
    assert!(centroid(&q).is_none());
    assert!(!intersects(&q, &rect(1, 1, 2, 2)));
    assert!(distance(&q, &p(50, 50)).is_none());
    assert!(!contains(&q, &rect(1, 1, 2, 2)));
}

/// Expected: `Shape::to_polygon` returns a counter-clockwise outer ring for any contour,
/// and picks the arc construction side from the true orientation.
/// Actual: for small contours away from the origin the outer ring comes back clockwise
/// (here a 3-unit triangle at 1 m = 1e9 nm; about half of all small random triangles near
/// `MAX_COORD` are affected), and arcs of such contours are approximated on the wrong side
/// (`Outside` <-> `Inside` swapped, since `material_left` is derived from the same sign).
///
/// Root cause: `contour_area` (src/arc.rs:231-247) sums `x0*y1 - x1*y0` in `f64` with
/// absolute coordinates: each term is ~1e18 (ulp 128) while the area is a few units, so
/// the sign is noise. The straight part could be summed exactly in `i128` relative to the
/// first vertex (like `ring_area2`), adding only the circular-segment terms in `f64`.
#[test]
fn shape_orientation_far_from_origin() {
    let pts = [
        p(1_000_767_959, 1_000_279_385),
        p(1_000_767_962, 1_000_279_385),
        p(1_000_767_958, 1_000_279_387),
    ];
    // Counter-clockwise as given (doubled area 6).
    assert_eq!(Ring::from(pts.to_vec()).signed_area2(), 6);
    let s = Shape::new(pts.iter().map(|&q| Curve::Line(q)).collect(), vec![]);
    let poly = s.to_polygon(ArcTol::new(10, Side::Outside)).unwrap();
    assert_eq!(poly.outer.signed_area2(), 6); // actual: -6 (reversed to clockwise)
}

/// Expected: the documented arc-reconstruction workflow (tag each arc, run booleans,
/// `arcs_from_tags` with the arc's centre and direction) gives back the result's
/// geometry. Actual: when the boolean traverses an arc backwards — a round pad subtracted
/// from a zone, i.e. every clearance cut-out, or any hole — the reconstructed
/// `CenterArc` uses the caller's original `ccw` and sweeps the complementary arc: here
/// the notch cut by the pad (traversed clockwise around the pad centre) is rebuilt as a
/// bulge on the other side of the chord (doubled area 1.95e10 becomes 2.01e10).
///
/// Root cause: src/arc.rs:614-631 (`arcs_from_tags`) emits `ccw` from `arc_of(tag)`
/// verbatim. The traversal direction is a property of the output ring, not of the tag;
/// it should be derived from the run itself (e.g. the sign of `orient(center, run start,
/// next vertex)`), or `arc_of` should return only the centre.
#[test]
fn arcs_from_tags_reversed_run() {
    let tol = ArcTol::new(10, Side::Outside);
    let pad = TaggedRing::uniform(Circle::new(p(0, 0), 10_000).to_ring(tol).unwrap(), 7);
    let zone = TaggedRing::uniform(rect(-5_000, -50_000, 100_000, 100_000), 0);
    let res = Boolean::new()
        .subject(&zone, FillRule::NonZero)
        .clip(&pad, FillRule::NonZero)
        .op(Op::Difference)
        .execute_tagged()
        .unwrap();
    assert_eq!(res.len(), 1);
    // The pad's arc as tagged: centre (0, 0), counter-clockwise (as approximated).
    let arc_of = |t: u64| (t == 7).then_some(p(0, 0));
    let contour = arcs_from_tags(&res[0].outer, &arc_of);
    let back = Shape::new(contour, vec![])
        .to_polygon(ArcTol::new(10, Side::Nearest))
        .unwrap();
    let orig = Ring(res[0].outer.points.clone()).signed_area2();
    let diff = (back.outer.signed_area2() - orig).abs();
    // Approximation differences are ~2 * tol * arc length ~ 1e6.
    assert!(
        diff < 10_000_000,
        "area2 {orig} rebuilt as {}",
        back.outer.signed_area2()
    );
}

/// Expected (documented in `triangulate`: "A vertex lying in the interior of another
/// ring's edge splits that edge" and "Two edges along the same segment in opposite
/// directions cancel each other (... polygons of a set sharing an edge)"): two rectangles
/// sharing part of an edge, or a hole sharing part of the outer ring's edge, triangulate.
/// Actual: `InvalidParameter("triangulate: invalid polygon ...")`. Sharing a whole edge
/// works; only the combination (split, then cancel) fails.
///
/// Root cause: opposite copies are cancelled once, before the sweep, on identical
/// `(lo, hi)` vertex pairs (src/triangulate.rs:272-290). Edges are split at vertices lying
/// on them only during the sweep (src/triangulate.rs:344-366), and the split piece is then
/// pushed next to its opposite copy, which the direction-order check
/// (src/triangulate.rs:374-378) rejects as overlapping.
#[test]
fn triangulate_partially_shared_edges() {
    let set = [
        Polygon::from(rect(0, 0, 20, 10)),
        Polygon::from(rect(0, 10, 10, 10)),
    ];
    assert_eq!(triangulate_set(&set).map(|t| t.area2()), Ok(2 * 300));
    let notch = Polygon::new(
        rect(0, 0, 10, 10),
        vec![Ring::from([(0, 2), (0, 6), (4, 6), (4, 2)])],
    );
    assert_eq!(triangulate(&notch).map(|t| t.area2()), Ok(2 * 84));
}

/// Expected: a boolean on `k` polygons that touch at one common vertex (e.g. re-unioning a
/// canonical result whose regions pinch at one point, which `trapezoids` and `offset` do
/// internally) runs in about `O(n log n)`. Actual: quadratic in the vertex degree: for
/// 8000 thin wedges around the origin (24 000 vertices) `union_all` of the canonical
/// union takes ~5.9 s (4000: 1.5 s, 16 000: ~20 s; `trapezoids` of 32 000 wedges ~68 s),
/// where the same number of disjoint triangles takes milliseconds.
///
/// Root cause (profiled): `pair_pass` / `for_each_pair` (src/node.rs:337-375): all `2k`
/// segments share the endpoint, so every leaf around it lists all of them and their
/// projections overlap in every direction: all O(k^2) pairs are tested. Then
/// `link_rings` (src/assemble.rs:88-103) scans every out-edge of a pinch vertex for each
/// arriving edge: O(k^2) again.
#[test]
#[cfg_attr(debug_assertions, ignore = "timing assertion: release builds only")]
fn union_of_polygons_sharing_one_vertex_is_fast() {
    let k = 8000;
    let r = 1_000_000_000.0f64;
    let wedges: Vec<Ring> = (0..k)
        .map(|i| {
            let a0 = std::f64::consts::PI * (2 * i) as f64 / k as f64;
            let a1 = std::f64::consts::PI * (2 * i + 1) as f64 / k as f64;
            Ring::from(vec![
                p(0, 0),
                p((r * a0.cos()) as i64, (r * a0.sin()) as i64),
                p((r * a1.cos()) as i64, (r * a1.sin()) as i64),
            ])
        })
        .collect();
    let u = union_all(&wedges, FillRule::NonZero).unwrap();
    assert_eq!(u.len(), k);
    let t = std::time::Instant::now();
    let again = union_all(&u, FillRule::NonZero).unwrap();
    let el = t.elapsed();
    assert_eq!(again, u);
    assert!(el.as_secs_f64() < 1.0, "re-union took {el:?}"); // actual ~5.9 s
}

/// Expected: `simplify_polygons` of a set does not depend on the order of its polygons
/// (the crate promises deterministic, canonical output; a set is unordered). Actual:
/// the result differs: `a`'s shortcut is blocked by `b`'s vertex (5, 10) when `a` is
/// processed first, but accepted when `b` was simplified first and that vertex is gone.
/// Ring rotation changes the result the same way (anchors are `pts[0]` and the vertex
/// farthest from it).
///
/// Root cause: src/simplify.rs:484-509 (`Simplifier::run`) runs Douglas-Peucker greedily
/// ring by ring in input order; accepted shortcuts change what later rings may do. A fix
/// would iterate to a fixpoint, or order the work canonically (e.g. sort rings by vertex
/// sequence and anchor at the minimum vertex).
#[test]
fn simplify_polygons_order_independent() {
    let a = Polygon::from(Ring::from([(0, 0), (10, 0), (10, 10), (5, 9), (0, 10)]));
    let b = Polygon::from(Ring::from([(0, 11), (5, 10), (10, 11), (10, 20), (0, 20)]));
    let ab = simplify_polygons(&[a.clone(), b.clone()], 2);
    let mut ba = simplify_polygons(&[b, a], 2);
    ba.reverse();
    // ab keeps (5, 9) in `a`; ba removes it.
    assert_eq!(ab, ba);
}

// =======================================================================================
// Regression coverage: investigated and found correct
// =======================================================================================

/// Several paths as one linear geometry (more than 512 segments in total, so `any_pair`
/// uses the direction-choosing sweep and `distance` the BVH).
struct Multi(Vec<Path>);

impl Geometry for Multi {
    fn bbox(&self) -> Option<Rect> {
        self.0
            .iter()
            .filter_map(|p| p.bbox())
            .reduce(|a, b| a.union(&b))
    }
    fn is_areal(&self) -> bool {
        false
    }
    fn visit_segments(&self, f: &mut dyn FnMut(Point, Point)) {
        for p in &self.0 {
            p.visit_segments(f)
        }
    }
    fn locate(&self, q: Point) -> Location {
        if self.0.iter().any(|p| p.locate(q) != Location::Outside) {
            Location::OnBoundary
        } else {
            Location::Outside
        }
    }
    fn any_point(&self) -> Option<Point> {
        self.0.first().and_then(|p| p.first().copied())
    }
}

/// O(n * m) reference: minimum over all segment pairs.
fn brute_distance(a: &[Path], b: &[Path]) -> SqDist {
    let segs = |x: &Path| -> Vec<Segment> {
        if x.len() == 1 {
            vec![Segment::new(x[0], x[0])]
        } else {
            x.windows(2).map(|w| Segment::new(w[0], w[1])).collect()
        }
    };
    let sb: Vec<Segment> = b.iter().flat_map(segs).collect();
    let mut best: Option<SqDist> = None;
    for u in a.iter().flat_map(segs) {
        for v in &sb {
            let s = distance_sq(&u, v).unwrap();
            if best.is_none_or(|b| s < b) {
                best = Some(s);
            }
        }
    }
    best.unwrap()
}

/// Dense parallel hatches at random angles (the case where the projection direction and
/// its margin matter), and random polylines spread over the whole coordinate range:
/// `distance`, `distance_less_than` and `intersects` agree with brute force.
#[test]
fn distance_queries_match_brute_force() {
    let mut rng = Rng(7);
    let m = MAX_COORD;
    for it in 0..24 {
        let (a, b) = if it % 2 == 0 {
            let hatch = |rng: &mut Rng, off: (i64, i64)| -> Vec<Path> {
                let (dx, dy) = (1 + rng.r(1000), rng.r(2001) - 1000);
                let (n, spacing, len) = (
                    150 + rng.r(250) as usize,
                    10 + rng.r(50),
                    3000 + rng.r(20000),
                );
                (0..n as i64)
                    .map(|i| {
                        let ox = off.0 - dy * i * spacing / 1000 + rng.r(3) - 1;
                        let oy = off.1 + dx * i * spacing / 1000 + rng.r(3) - 1;
                        Path::from(vec![
                            p(ox, oy),
                            p(ox + dx * len / 1000, oy + dy * len / 1000),
                        ])
                    })
                    .collect()
            };
            let a = hatch(&mut rng, (0, 0));
            let off = (rng.r(20000) - 10000, rng.r(20000) - 10000);
            (a, hatch(&mut rng, off))
        } else {
            let scale = [m, m / 1000, 1 << 20][(it / 2) % 3];
            let mk = |rng: &mut Rng, ox: i64| -> Vec<Path> {
                let n = 150 + rng.r(200) as usize;
                (0..n)
                    .map(|_| {
                        let (x, y) = (rng.r(2 * scale) - scale, rng.r(2 * scale) - scale);
                        let mut v = vec![p((x / 2 + ox).clamp(-m, m), y)];
                        for _ in 0..1 + rng.r(4) {
                            let l = *v.last().unwrap();
                            let dx = rng.r(scale / 8) - scale / 16;
                            let dy = rng.r(scale / 8) - scale / 16;
                            v.push(p((l.x + dx).clamp(-m, m), (l.y + dy).clamp(-m, m)));
                        }
                        Path(v)
                    })
                    .collect()
            };
            let a = mk(&mut rng, -scale / 2);
            (a, mk(&mut rng, scale / 2))
        };
        let bf = brute_distance(&a, &b);
        let (a, b) = (Multi(a), Multi(b));
        assert_eq!(distance(&a, &b).unwrap().sq, bf, "it {it}");
        assert_eq!(intersects(&a, &b), bf.is_zero(), "it {it}");
        let df = bf.distance_f64();
        for d in [df.floor() as i64, df.ceil() as i64, df.ceil() as i64 + 1] {
            let exp = d > 0 && bf.cmp_dist(d as u64) == Ordering::Less;
            assert_eq!(distance_less_than(&a, &b, d), exp, "it {it}: d {d}");
        }
    }
}

fn random_set(rng: &mut Rng, range: i64) -> PolygonSet {
    let rings: Vec<Ring> = (0..1 + rng.r(8))
        .map(|_| {
            if rng.r(2) == 0 {
                let (x, y) = (rng.r(range), rng.r(range));
                rect(x, y, 1 + rng.r(range / 2), 1 + rng.r(range / 2))
            } else {
                (0..3 + rng.r(8))
                    .map(|_| p(rng.r(range), rng.r(range)))
                    .collect()
            }
        })
        .collect();
    let rule = [FillRule::NonZero, FillRule::EvenOdd][rng.r(2) as usize];
    union_all(&rings, rule).unwrap()
}

/// `trapezoids` covers the region exactly with disjoint interiors: total area matches, no
/// sample point (on a finer grid) lies strictly inside two trapezoids or inside a
/// trapezoid but outside the region, and every trapezoid's edges span `[x0, x1]`.
#[test]
fn trapezoids_partition_the_region() {
    let mut rng = Rng(11);
    let orient = |a: Point, b: Point, c: Point| -> i128 {
        (b.x - a.x) as i128 * (c.y - a.y) as i128 - (b.y - a.y) as i128 * (c.x - a.x) as i128
    };
    for it in 0..600 {
        let range = [10, 30, 200][it % 3];
        let set = random_set(&mut rng, range);
        let tr = trapezoids(&set).unwrap();
        let area: f64 = tr.iter().map(|t| t.area()).sum();
        let exp = area2(&set) as f64 / 2.0;
        assert!((area - exp).abs() < 1e-6 * exp.max(1.0), "it {it}");
        for t in &tr {
            assert!(t.x1 > t.x0, "it {it}");
            assert!(t.bottom.0.x <= t.x0 && t.bottom.1.x >= t.x1, "it {it}");
            assert!(t.top.0.x <= t.x0 && t.top.1.x >= t.x1, "it {it}");
        }
        // Sample on a grid 7 times finer (scaled coordinates).
        let s = |a: Point| p(a.x * 7, a.y * 7);
        let scaled: PolygonSet = set
            .iter()
            .map(|q| Polygon {
                outer: q.outer.iter().map(|&v| s(v)).collect(),
                holes: q
                    .holes
                    .iter()
                    .map(|h| h.iter().map(|&v| s(v)).collect())
                    .collect(),
            })
            .collect();
        for _ in 0..100 {
            let q = p(rng.r(range * 21 / 2), rng.r(range * 21 / 2));
            let c = tr
                .iter()
                .filter(|t| {
                    q.x > t.x0 * 7
                        && q.x < t.x1 * 7
                        && orient(s(t.bottom.0), s(t.bottom.1), q) > 0
                        && orient(s(t.top.0), s(t.top.1), q) < 0
                })
                .count();
            if locate(&scaled, q) == Location::Inside {
                assert!(c <= 1, "it {it}: {q:?} in {c} trapezoids");
            } else {
                assert_eq!(c, 0, "it {it}: {q:?} outside the region");
            }
        }
    }
}

/// Triangulations (plain and Delaunay) of random canonical polygons (with holes and pinch
/// vertices) cover the polygon area, do not depend on ring orientation, start vertex or
/// hole order, and the Delaunay result is locally Delaunay on every interior edge.
#[test]
fn triangulations_are_invariant_and_delaunay() {
    fn incircle(a: Point, b: Point, c: Point, d: Point) -> i128 {
        let f = |p: Point| ((p.x - d.x) as i128, (p.y - d.y) as i128);
        let ((ax, ay), (bx, by), (cx, cy)) = (f(a), f(b), f(c));
        let (a2, b2, c2) = (ax * ax + ay * ay, bx * bx + by * by, cx * cx + cy * cy);
        ax * (by * c2 - b2 * cy) - ay * (bx * c2 - b2 * cx) + a2 * (bx * cy - by * cx)
    }
    let mut rng = Rng(12);
    for it in 0..600 {
        let range = [10, 30, 200][it % 3];
        for poly in &random_set(&mut rng, range) {
            let t = triangulate(poly).unwrap();
            let d = triangulate_delaunay(poly).unwrap();
            assert_eq!(t.area2(), poly.signed_area2(), "it {it}");
            assert_eq!(d.area2(), poly.signed_area2(), "it {it}");
            let mut q = poly.clone();
            let l = q.outer.len();
            q.outer.rotate_left(it % l);
            q.outer.reverse();
            q.holes.iter_mut().for_each(|h| h.reverse());
            q.holes.reverse();
            assert_eq!(triangulate(&q).unwrap(), t, "it {it}");
            assert_eq!(triangulate_delaunay(&q).unwrap(), d, "it {it}");
            let boundary: std::collections::HashSet<(Point, Point)> = poly
                .rings()
                .flat_map(|r| r.edges())
                .flat_map(|(a, b)| [(a, b), (b, a)])
                .collect();
            let tris: Vec<[Point; 3]> = (0..d.triangles.len()).map(|i| d.triangle(i)).collect();
            let mut by_edge = std::collections::HashMap::new();
            for (i, tr) in tris.iter().enumerate() {
                for k in 0..3 {
                    by_edge.insert((tr[k], tr[(k + 1) % 3]), i);
                }
            }
            for tr in &tris {
                for k in 0..3 {
                    let (a, b, c) = (tr[k], tr[(k + 1) % 3], tr[(k + 2) % 3]);
                    if boundary.contains(&(a, b)) {
                        continue;
                    }
                    if let Some(&j) = by_edge.get(&(b, a)) {
                        let o = *tris[j].iter().find(|&&v| v != a && v != b).unwrap();
                        assert!(incircle(a, b, c, o) <= 0, "it {it}: {a:?}-{b:?}");
                    }
                }
            }
        }
    }
}

/// Round-join offsets agree with exact distances on the requested side: for `Outside`
/// every point within `|delta| - 1` of the region is covered (growing) or every point
/// deeper than `|delta| + 1` is kept (shrinking), and symmetrically for `Inside`; small
/// deltas exercise the concave-join trimming.
#[test]
fn round_offsets_respect_the_side() {
    let mut rng = Rng(21);
    for it in 0..200 {
        let set = random_set(&mut rng, 10_000);
        let d = if it % 2 == 0 {
            1 + rng.r(2000)
        } else {
            1 + rng.r(20)
        };
        let delta = if rng.r(2) == 0 { d } else { -d };
        let tol = 1 + rng.r(6);
        let side = [Side::Outside, Side::Inside, Side::Nearest][rng.r(3) as usize];
        let off = offset(&set, delta, Join::Round, ArcTol::new(tol, side)).unwrap();
        assert_eq!(check_canonical(&off, true), Ok(()), "it {it}");
        let t = tol as f64 + 1.0;
        let (bin, bout) = match side {
            Side::Outside => (1.0, t),
            Side::Inside => (t, 1.0),
            Side::Nearest => (t, t),
        };
        let boundary: Vec<Path> = set
            .iter()
            .flat_map(|q| q.rings())
            .map(|r| r.iter().copied().chain([r[0]]).collect())
            .collect();
        for _ in 0..60 {
            let q = p(rng.r(15_000) - 2_500, rng.r(15_000) - 2_500);
            let inside = locate(&set, q) != Location::Outside;
            let db = boundary
                .iter()
                .map(|b| distance_sq(b, &q).unwrap().distance_f64())
                .fold(f64::INFINITY, f64::min);
            let ad = d as f64;
            let expect = if delta > 0 {
                if inside || db < ad - bin {
                    Some(true)
                } else if db > ad + bout {
                    Some(false)
                } else {
                    None
                }
            } else if !inside || db < ad - bout {
                Some(false)
            } else if db > ad + bin {
                Some(true)
            } else {
                None
            };
            if let Some(e) = expect {
                assert_eq!(locate(&off, q) != Location::Outside, e, "it {it}: {q:?}");
            }
        }
    }
}
