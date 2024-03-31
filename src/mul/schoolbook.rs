use super::*;
use equator::debug_assert;

fn widening_mul_acc(lhs: Limb, rhs: Limb, acc: Limb, carry: Limb) -> (Limb, Limb) {
    let wide = acc as u128 + carry as u128 + lhs as u128 * rhs as u128;
    (wide as u64, (wide >> 64) as u64)
}

#[inline]
pub fn mul_bigint(full_mul: &mut [Limb], lhs: &[Limb], rhs: &[Limb]) {
    debug_assert!(full_mul.len() == lhs.len() + rhs.len());
    full_mul.fill(consts::LIMB_ZERO);
    for (i, &l) in lhs.iter().enumerate() {
        let dst = &mut full_mul[i..];
        let mut carry = consts::LIMB_ZERO;
        for (dst, &r) in core::iter::zip(&mut *dst, rhs) {
            (*dst, carry) = widening_mul_acc(l, r, *dst, carry);
        }
        dst[rhs.len()] = carry;
    }
}
