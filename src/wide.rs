//! Minimal unsigned wide-integer arithmetic for exact comparisons of large products.

use core::cmp::Ordering;

/// Unsigned 384-bit integer, little-endian 64-bit limbs.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub(crate) struct U384(pub [u64; 6]);

impl U384 {
    #[inline]
    pub fn from_u128(v: u128) -> Self {
        U384([v as u64, (v >> 64) as u64, 0, 0, 0, 0])
    }

    /// Product of two `u128` (fits in 256 bits).
    #[inline]
    pub fn mul_u128(a: u128, b: u128) -> Self {
        if (a >> 64) == 0 && (b >> 64) == 0 {
            return Self::from_u128(a * b);
        }
        Self::from_u128(a).mul_small(b)
    }

    /// `self * b`, truncated to 384 bits (callers keep products below 2^384).
    pub fn mul_small(&self, b: u128) -> Self {
        let bl = [b as u64, (b >> 64) as u64];
        let mut r = [0u64; 6];
        for (j, &bj) in bl.iter().enumerate() {
            if bj == 0 {
                continue;
            }
            let mut carry: u128 = 0;
            for i in 0..6 {
                if i + j >= 6 {
                    break;
                }
                let t = self.0[i] as u128 * bj as u128 + r[i + j] as u128 + carry;
                r[i + j] = t as u64;
                carry = t >> 64;
            }
        }
        U384(r)
    }

    /// Converts to `f64` (rounded).
    pub fn to_f64(self) -> f64 {
        let mut v = 0f64;
        for i in (0..6).rev() {
            v = v * 18446744073709551616.0 + self.0[i] as f64;
        }
        v
    }
}

impl Ord for U384 {
    fn cmp(&self, o: &Self) -> Ordering {
        for i in (0..6).rev() {
            match self.0[i].cmp(&o.0[i]) {
                Ordering::Equal => {}
                x => return x,
            }
        }
        Ordering::Equal
    }
}

impl PartialOrd for U384 {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}

/// Compares `a * b` with `c * d` exactly.
#[inline]
pub(crate) fn cmp_products(a: u128, b: u128, c: u128, d: u128) -> Ordering {
    if let (Some(x), Some(y)) = (a.checked_mul(b), c.checked_mul(d)) {
        return x.cmp(&y);
    }
    U384::mul_u128(a, b).cmp(&U384::mul_u128(c, d))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn products() {
        let a = u128::MAX;
        let p = U384::mul_u128(a, a);
        // (2^128 - 1)^2 = 2^256 - 2^129 + 1
        assert_eq!(p.0, [1, 0, u64::MAX - 1, u64::MAX, 0, 0]);
        assert_eq!(cmp_products(a, 2, 2, a), Ordering::Equal);
        assert_eq!(cmp_products(a, 3, 2, a), Ordering::Greater);
        let q = p.mul_small(1 << 100);
        assert!(q > p);
        assert_eq!(U384::from_u128(12345).to_f64(), 12345.0);
    }
}
