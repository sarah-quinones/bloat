use super::*;
use equator::debug_assert;

#[inline]
pub fn mul_bigint(full_mul: &mut [Limb], lhs: &[Limb], rhs: &[Limb]) {
    debug_assert!(full_mul.len() == lhs.len() + rhs.len());
    full_mul.fill(consts::LIMB_ZERO);
    for (i, &l) in lhs.iter().enumerate() {
        let dst = &mut full_mul[i..];
        let mut carry = consts::LIMB_ZERO;
        for (dst, &r) in core::iter::zip(&mut *dst, rhs) {
            let big = *dst as u128 + carry as u128 + r as u128 * l as u128;
            *dst = big as u64;
            carry = (big >> 64) as u64;
        }
        dst[rhs.len()] = carry;
    }
}
