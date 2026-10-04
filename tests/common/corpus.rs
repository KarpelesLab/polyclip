//! Reader for the real-world corpus in `testdata/` (`.pclp` files): polygon sets dumped
//! from cadlab boards, stored compactly as zigzag-varint coordinate deltas.
//!
//! Format: `b"PCLP1\n"`, then varints: polygon count; per polygon its ring count (outer
//! first, then holes); per ring its vertex count and the vertices as deltas from the
//! previous vertex of the same ring (starting from the origin).

#![allow(dead_code)]

use polyclip::{Point, Polygon, PolygonSet, Ring};

/// Loads a corpus file from `testdata/` (`name` without directory).
pub fn load(name: &str) -> PolygonSet {
    let path = format!("{}/testdata/{name}", env!("CARGO_MANIFEST_DIR"));
    let data = std::fs::read(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    assert!(data.starts_with(b"PCLP1\n"), "{path}: bad header");
    let mut pos = 6;
    let mut next = || -> i64 {
        let (mut v, mut shift) = (0u64, 0);
        loop {
            let b = data[pos];
            pos += 1;
            v |= ((b & 0x7f) as u64) << shift;
            if b & 0x80 == 0 {
                break;
            }
            shift += 7;
        }
        ((v >> 1) as i64) ^ -((v & 1) as i64)
    };
    let np = next() as usize;
    (0..np)
        .map(|_| {
            let nr = next() as usize;
            let mut rings: Vec<Ring> = (0..nr)
                .map(|_| {
                    let n = next() as usize;
                    let (mut x, mut y) = (0i64, 0i64);
                    (0..n)
                        .map(|_| {
                            x += next();
                            y += next();
                            Point::new(x, y)
                        })
                        .collect()
                })
                .collect();
            let outer = rings.remove(0);
            Polygon::new(outer, rings)
        })
        .collect()
}
