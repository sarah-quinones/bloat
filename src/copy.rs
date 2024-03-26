use super::*;
use equator::assert;

fn copy_finite_with_sign_round<L: Limb>(
    dst: &mut BigFloat<L>,
    src: &BigFloat<L>,
    rnd: Round,
    sign: Sign,
) -> Approx {
    assert!(src.precision_bits() > dst.precision_bits());
    assert!(src.mantissa().len() >= dst.mantissa().len());

    dst.sign_biased_exponent = ((sign.is_negative() as u64) << consts::SIGN_SHIFT)
        | (src.sign_biased_exponent & consts::BIASED_EXPONENT_MASK);

    match src.exponent() {
        Exponent::Zero | Exponent::NaN | Exponent::Inf => return Approx::Exact,
        _ => {}
    }

    let limb_diff = src.mantissa().len() - dst.mantissa().len();

    let rnd = rnd.with_sign(sign);

    let one = L::ONE;

    let dst_bits_mod_limb = dst.precision_bits() % L::BITS;
    let amount = L::BITS - dst_bits_mod_limb - 1;

    let (msb, lsb_any) = if dst_bits_mod_limb == 0 {
        let msb = src.mantissa()[limb_diff - 1].shr(amount) == one;
        let mut lsb_any = src.mantissa()[limb_diff - 1] & (one.shl(amount).wrapping_sub(one));
        for &l in src.mantissa()[..limb_diff - 1].iter().rev() {
            lsb_any |= l;
        }
        (msb, lsb_any != L::ZERO)
    } else {
        let msb = (src.mantissa()[limb_diff].shr(amount) & one) == one;
        let mut lsb_any = src.mantissa()[limb_diff] & (one.shl(amount).wrapping_sub(one));
        for &l in src.mantissa()[..limb_diff].iter().rev() {
            lsb_any |= l;
        }
        (msb, lsb_any != L::ZERO)
    };

    if !msb && !lsb_any {
        dst.mantissa_mut()
            .copy_from_slice(&src.mantissa()[limb_diff..]);
        Approx::Exact
    } else {
        let amount = (L::BITS - dst_bits_mod_limb) % L::BITS;
        let ulp = one.shl(amount);
        let mantissa_odd = (src.mantissa()[limb_diff].shr(amount) & one) == one;

        if rnd == RoundKnownSign::AwayFromZero
            || (rnd == RoundKnownSign::ToNearest && msb && (lsb_any || mantissa_odd))
        {
            // round away from zero
            let mut carry;
            let src_m = &src.mantissa()[limb_diff..];

            (dst.mantissa_mut()[0], carry) =
                (src_m[0] & !(ulp.wrapping_sub(one))).overflowing_add(ulp);

            for (dst, &src) in core::iter::zip(&mut dst.mantissa_mut()[1..], &src_m[1..]) {
                (*dst, carry) = src.overflowing_add(L::from_bit(carry));
            }

            if carry {
                dst.sign_biased_exponent += 1;
            }

            Approx::from_sign(sign)
        } else {
            // round to zero
            dst.mantissa_mut()
                .copy_from_slice(&src.mantissa()[limb_diff..]);
            dst.mantissa_mut()[0] &= !(ulp.wrapping_sub(one));

            Approx::from_sign(sign.neg())
        }
    }
}

pub fn copy_with_sign<L: Limb>(
    dst: &mut BigFloat<L>,
    src: &BigFloat<L>,
    rnd: Round,
    sign: Sign,
) -> Approx {
    dst.sign_biased_exponent = ((sign.is_negative() as u64) << consts::SIGN_SHIFT)
        | (src.sign_biased_exponent & consts::BIASED_EXPONENT_MASK);
    let approx = if dst.precision_bits() >= src.precision_bits() {
        dst.mantissa_mut()[..src.mantissa().len()].copy_from_slice(&src.mantissa());
        dst.mantissa_mut()[src.mantissa().len()..].fill(L::ZERO);
        Approx::Exact
    } else {
        copy_finite_with_sign_round(dst, src, rnd, sign)
    };
    approx
}

#[inline]
pub fn copy<L: Limb>(dst: &mut BigFloat<L>, src: &BigFloat<L>, rnd: Round) -> Approx {
    copy_with_sign(dst, src, rnd, src.sign())
}

#[inline]
pub fn abs<L: Limb>(dst: &mut BigFloat<L>, src: &BigFloat<L>, rnd: Round) -> Approx {
    copy_with_sign(dst, src, rnd, Sign::Pos)
}

#[cfg(test)]
mod tests {
    use super::*;
    use equator::assert;

    #[test]
    fn test_widening_copy() {
        let src = SmallFloat::from_parts(
            Sign::Neg,
            Exponent::Finite(3),
            4,
            utils::rev([0b1011_0000_0000_0000_0000_0000_0000_0000_u32]),
        );
        let mut dst = SmallFloat::from_parts(Sign::Pos, Exponent::Zero, 8, utils::rev([0_u32]));

        for &rnd in Round::ALL {
            assert!(copy(&mut dst, &src, rnd) == Approx::Exact);

            assert!(all(
                dst.precision_bits() == 8,
                dst.exponent() == Exponent::Finite(3),
                dst.sign() == Sign::Neg,
                dst.mantissa() == &[0b1011_0000_0000_0000_0000_0000_0000_0000_u32]
            ));
        }
    }

    #[test]
    fn test_very_widening_copy() {
        let src = SmallFloat::from_parts(
            Sign::Neg,
            Exponent::Finite(3),
            4,
            utils::rev([0b1011_0000_0000_0000_0000_0000_0000_0000_u32]),
        );
        let mut dst = SmallFloat::from_parts(Sign::Pos, Exponent::Zero, 40, utils::rev([0, 0_u32]));

        for &rnd in Round::ALL {
            assert!(copy(&mut dst, &src, rnd) == Approx::Exact);

            assert!(all(
                dst.precision_bits() == 40,
                dst.exponent() == Exponent::Finite(3),
                dst.sign() == Sign::Neg,
                dst.mantissa()
                    == &[
                        0b1011_0000_0000_0000_0000_0000_0000_0000,
                        0b0000_0000_0000_0000_0000_0000_0000_0000_u32
                    ]
            ));
        }
    }

    #[test]
    fn test_narrowing_copy_exact() {
        let src = SmallFloat::from_parts(
            Sign::Neg,
            Exponent::Finite(3),
            8,
            utils::rev([0b1011_0000_0000_0000_0000_0000_0000_0000_u32]),
        );
        let mut dst = SmallFloat::from_parts(Sign::Pos, Exponent::Zero, 4, utils::rev([0, 0_u32]));

        for &rnd in Round::ALL {
            assert!(copy(&mut dst, &src, rnd) == Approx::Exact);

            assert!(all(
                dst.precision_bits() == 4,
                dst.exponent() == Exponent::Finite(3),
                dst.sign() == Sign::Neg,
                dst.mantissa() == &[0b1011_0000_0000_0000_0000_0000_0000_0000]
            ));
        }
    }

    #[test]
    fn test_narrowing_copy_round() {
        let src = SmallFloat::from_parts(
            Sign::Neg,
            Exponent::Finite(3),
            8,
            utils::rev([0b1011_0000_0000_0000_0000_0000_0000_0000_u32]),
        );
        let mut dst = SmallFloat::from_parts(Sign::Pos, Exponent::Zero, 4, utils::rev([0, 0_u32]));

        for &rnd in Round::ALL {
            assert!(copy(&mut dst, &src, rnd) == Approx::Exact);
            assert!(all(
                dst.precision_bits() == 4,
                dst.exponent() == Exponent::Finite(3),
                dst.sign() == Sign::Neg,
                dst.mantissa() == &[0b1011_0000_0000_0000_0000_0000_0000_0000]
            ));
        }
    }
}
