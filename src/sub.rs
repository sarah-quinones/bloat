use super::*;

fn sign_exponent_of_add_different_sign(lhs_sign: Sign, lhs: &BigFloat, rhs: &BigFloat) -> (Sign, Exponent, bool) {
    let mut underflow = false;

    let (mut swap, lhs, rhs) =
        if (lhs.sign_biased_exponent & consts::BIASED_EXPONENT_MASK) >= (rhs.sign_biased_exponent & consts::BIASED_EXPONENT_MASK) {
            (false, lhs, rhs)
        } else {
            (true, rhs, lhs)
        };

    let lhs_exp = lhs.exponent();
    let rhs_exp = rhs.exponent();
    if rhs_exp == Exponent::Inf || rhs_exp == Exponent::NaN || lhs_exp == Exponent::NaN {
        return (Sign::Pos, Exponent::NaN, underflow);
    }
    if lhs_exp == Exponent::Inf || lhs_exp == Exponent::Zero || rhs_exp == Exponent::Zero {
        let sign = lhs_sign;
        return (if swap { sign.neg() } else { sign }, lhs_exp, underflow);
    }

    let (Exponent::Finite(lhs_exp), Exponent::Finite(rhs_exp)) = (lhs_exp, rhs_exp) else {
        unreachable!()
    };

    let mut lhs = lhs.mantissa();
    let mut rhs = rhs.mantissa();

    let mut res = consts::LIMB_ZERO;

    // remove the common limbs, if any
    if lhs_exp == rhs_exp {
        let mut equal = true;

        while let (Some((&l, lhs_head)), Some((&r, rhs_head))) = (lhs.split_last(), rhs.split_last()) {
            if l != r {
                if l < r {
                    swap = !swap;
                    (lhs, rhs) = (rhs, lhs);
                }
                equal = false;
                break;
            }
            res += consts::LIMB_BITS;
            lhs = lhs_head;
            rhs = rhs_head;
        }

        if equal {
            if lhs.len() == 0 {
                swap = !swap;
                (lhs, _) = (rhs, lhs);
            }

            for &l in lhs.iter().rev() {
                if l != consts::LIMB_ZERO {
                    res += l.leading_zeros() as u64;

                    let exp = lhs_exp.saturating_sub_unsigned(res);
                    let sign = if swap { Sign::Neg } else { Sign::Pos };
                    let exp = if exp < consts::MIN_EXPONENT_INCLUSIVE {
                        underflow = true;
                        Exponent::Zero
                    } else {
                        Exponent::Finite(exp)
                    };
                    return (sign, exp, underflow);
                }
                res += consts::LIMB_BITS;
            }

            return (Sign::Pos, Exponent::Zero, underflow);
        }
    }

    let exp_shift = (lhs_exp as u64).wrapping_sub(rhs_exp as u64);
    let exp_large_shift = (exp_shift / consts::LIMB_BITS) as usize;
    let exp_small_shift = exp_shift % consts::LIMB_BITS;

    let zero = &consts::LIMB_ZERO;

    let lhs_len = lhs.len();
    let rhs_len = rhs.len() + exp_large_shift;
    let lhs = |i: usize| {
        *if i < lhs.len() { &lhs[lhs.len() - i - 1] } else { zero }
    };

    let rhs_large_shift = |i: usize| {
        *if i.wrapping_sub(exp_large_shift) < rhs.len() {
            &rhs[exp_large_shift + rhs.len() - i - 1]
        } else {
            zero
        }
    };
    let rhs = |i: usize| {
        let r0 = rhs_large_shift(i);

        if exp_small_shift == 0 || i == 0 {
            r0.shr(exp_small_shift)
        } else {
            let r1 = rhs_large_shift(i - 1);
            r0.shr(exp_small_shift) | r1.shl(consts::LIMB_BITS - exp_small_shift)
        }
    };

    let mut i = 0;

    let mut dif_hi = consts::LIMB_ZERO;
    let mut dif_lo = lhs(i).wrapping_sub(rhs(i));

    while dif_lo == consts::LIMB_ONE && dif_hi == consts::LIMB_ZERO {
        i += 1;
        let (new_dif_lo, borrow) = lhs(i).overflowing_sub(rhs(i));
        dif_hi = dif_lo.wrapping_sub(borrow as Limb);
        dif_lo = new_dif_lo;
        res += consts::LIMB_BITS;
    }

    if dif_hi == consts::LIMB_ONE {
        res -= 1;
    } else {
        res += dif_lo.leading_zeros() as u64;
    }

    let pow_of_2 =
        (dif_hi == consts::LIMB_ONE && dif_lo == consts::LIMB_ZERO) || (dif_lo & dif_lo.wrapping_sub(consts::LIMB_ONE) == consts::LIMB_ZERO);

    if pow_of_2 {
        loop {
            i += 1;
            if i >= Ord::max(lhs_len, rhs_len) || lhs(i) != rhs(i) {
                break;
            }
        }
        res += (lhs(i) < rhs(i)) as u64;
    }

    let exp = lhs_exp.saturating_sub_unsigned(res);
    let sign = if swap { lhs_sign.neg() } else { lhs_sign };
    let exp = if exp < consts::MIN_EXPONENT_INCLUSIVE {
        underflow = true;
        Exponent::Zero
    } else {
        Exponent::Finite(exp)
    };
    (sign, exp, underflow)
}

pub fn add_different_sign(dst: &mut BigFloat, lhs_sign: Sign, lhs: &BigFloat, rhs: &BigFloat, rnd: Round) -> Approx {
    let mut lhs = lhs;
    let mut rhs = rhs;

    let (sign, dst_exp, underflow) = sign_exponent_of_add_different_sign(lhs_sign, lhs, rhs);
    dst.sign_biased_exponent = crate::make_sign_and_biased_exponent(sign, dst_exp);
    let dst_exp = match dst_exp {
        Exponent::Zero if underflow => return Approx::Underflow,
        Exponent::Zero | Exponent::Inf | Exponent::NaN => return Approx::Exact,
        Exponent::Finite(exp) => exp,
    };

    if sign.is_negative() {
        (lhs, rhs) = (rhs, lhs)
    }

    if rhs.exponent() == Exponent::Zero {
        dst.copy_from(lhs, rnd);
    }

    let (Exponent::Finite(lhs_exp), Exponent::Finite(rhs_exp)) = (lhs.exponent(), rhs.exponent()) else {
        unreachable!()
    };

    let rnd = rnd.with_sign(sign);

    let lhs_exp_shift = (lhs_exp as u64).wrapping_sub(dst_exp as u64);
    let rhs_exp_shift = rhs_exp - dst_exp; // may be negative

    let lhs_exp_large_shift = (lhs_exp_shift / consts::LIMB_BITS) as usize;
    let lhs_exp_small_shift = lhs_exp_shift % consts::LIMB_BITS;

    let rhs_exp_small_shift = rhs_exp_shift.rem_euclid(consts::LIMB_BITS as i64);
    let rhs_exp_large_shift = ((rhs_exp_shift - rhs_exp_small_shift) / consts::LIMB_BITS as i64) as isize;
    let rhs_exp_small_shift = rhs_exp_small_shift as u64;

    let lhs = lhs.mantissa();
    let rhs = rhs.mantissa();

    let lhs_len = lhs.len();
    let rhs_len = rhs.len();

    let zero = &consts::LIMB_ZERO;

    let lhs_large_shift = |i: usize| {
        let idx = i.wrapping_add(lhs_exp_large_shift);
        *if idx < lhs.len() { &lhs[lhs.len() - idx - 1] } else { zero }
    };
    let rhs_large_shift = |i: usize| {
        let idx = i.wrapping_add(rhs_exp_large_shift as usize);
        *if idx < rhs.len() { &rhs[rhs.len() - idx - 1] } else { zero }
    };

    let lhs = |i: usize| {
        let l0 = lhs_large_shift(i);
        let l1 = lhs_large_shift(i + 1);

        if lhs_exp_small_shift == 0 {
            l0
        } else {
            l0.shl(lhs_exp_small_shift) | l1.shr(consts::LIMB_BITS - lhs_exp_small_shift)
        }
    };
    let rhs = |i: usize| {
        let r0 = rhs_large_shift(i);
        let r1 = rhs_large_shift(i + 1);

        if rhs_exp_small_shift == 0 {
            r0
        } else {
            r0.shl(rhs_exp_small_shift) | r1.shr(consts::LIMB_BITS - rhs_exp_small_shift)
        }
    };

    let one = consts::LIMB_ONE;

    let dst_len = dst.mantissa_len();
    let dst_bits_mod_limb = dst.precision_bits() % consts::LIMB_BITS;

    let idx_that_contains_msb = if dst_bits_mod_limb == 0 { dst_len } else { dst_len - 1 };
    let msb_pos = consts::LIMB_BITS - dst_bits_mod_limb - 1;
    let ulp_pos = (consts::LIMB_BITS - dst_bits_mod_limb) % consts::LIMB_BITS;
    let ulp = one.shl(ulp_pos);

    let dst_m = dst.mantissa_mut();
    let max = Ord::max;
    let n = max(
        dst_len,
        max(rhs_len.saturating_add_signed(rhs_exp_large_shift), lhs_len + lhs_exp_large_shift),
    );

    let mut lsb_any_word = consts::LIMB_ZERO;

    let mut borrow = false;
    for i in (idx_that_contains_msb + 1..n).rev() {
        let (diff, borrow0) = lhs(i).overflowing_sub(rhs(i));
        let (diff, borrow1) = diff.overflowing_sub(borrow as Limb);
        borrow = borrow0 | borrow1;

        lsb_any_word |= diff;
    }
    let i = idx_that_contains_msb;
    let (diff, borrow0) = lhs(i).overflowing_sub(rhs(i));
    let (diff, borrow1) = diff.overflowing_sub(borrow as Limb);
    borrow = borrow0 | borrow1;

    lsb_any_word |= diff & one.shl(msb_pos).wrapping_sub(one);
    let msb = (diff & one.shl(msb_pos)) == one.shl(msb_pos);
    let lsb_any = lsb_any_word != consts::LIMB_ZERO;

    if i == dst_len - 1 {
        dst_m[0] = diff;
    }

    let result_lsb = (dst_m[0] & ulp) == ulp;

    for i in (0..idx_that_contains_msb).rev() {
        let (diff, borrow0) = lhs(i).overflowing_sub(rhs(i));
        let (diff, borrow1) = diff.overflowing_sub(borrow as Limb);
        borrow = borrow0 | borrow1;
        let idx = dst_len - i - 1;
        dst_m[idx] = diff;
    }

    if rnd == RoundKnownSign::AwayFromZero || (rnd == RoundKnownSign::ToNearest && msb && (lsb_any || result_lsb)) {
        // add 1ulp
        let mut carry;
        (dst_m[0], carry) = dst_m[0].overflowing_add(ulp);
        for x in &mut dst_m[1..] {
            (*x, carry) = (*x).overflowing_add(carry as Limb);
        }

        if carry {
            dst.sign_biased_exponent += 1;
            *dst.mantissa_mut().last_mut().unwrap() = one.shl(consts::LIMB_BITS - 1);
        }
    }
    dst.mantissa_mut()[0] &= !(one.shl(ulp_pos).wrapping_sub(one));

    if !msb && !lsb_any {
        Approx::Exact
    } else if rnd == RoundKnownSign::AwayFromZero || (rnd == RoundKnownSign::ToNearest && msb && (lsb_any || result_lsb)) {
        if dst.sign_biased_exponent & consts::BIASED_EXPONENT_MASK == consts::BIASED_EXPONENT_INF {
            Approx::Overflow
        } else {
            Approx::from_sign(sign)
        }
    } else {
        Approx::from_sign(sign.neg())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use equator::assert;

    #[test]
    fn test_sign_exponent_of_add_0() {
        let a = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(2),
            utils::rev([0b1011_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let b = SmallFloat::from_parts(
            4,
            Sign::Neg,
            Exponent::Finite(2),
            utils::rev([0b1001_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        assert!(sign_exponent_of_add_different_sign(a.sign(), &a, &b) == (Sign::Pos, Exponent::Finite(0), false));
    }

    #[test]
    fn test_sign_exponent_of_add_1() {
        let a = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(2),
            utils::rev([0b1011_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let b = SmallFloat::from_parts(
            4,
            Sign::Neg,
            Exponent::Finite(2),
            utils::rev([0b1010_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        assert!(sign_exponent_of_add_different_sign(a.sign(), &a, &b) == (Sign::Pos, Exponent::Finite(-1), false));
    }

    #[test]
    fn test_sign_exponent_of_add_2() {
        let a = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(2),
            utils::rev([0b1010_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let b = SmallFloat::from_parts(
            4,
            Sign::Neg,
            Exponent::Finite(2),
            utils::rev([0b1111_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        assert!(sign_exponent_of_add_different_sign(a.sign(), &a, &b) == (Sign::Neg, Exponent::Finite(1), false));
    }

    #[test]
    fn test_sub_0() {
        let a = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(2),
            utils::rev([0b1011_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let b = SmallFloat::from_parts(
            4,
            Sign::Neg,
            Exponent::Finite(2),
            utils::rev([0b1001_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let mut c = SmallFloat::<1>::zero(4);
        add_different_sign(&mut c, a.sign(), &a, &b, Round::ToNearest);
        assert!(
            c.repr()
                == SmallFloat::from_parts(
                    4,
                    Sign::Pos,
                    Exponent::Finite(0),
                    utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
                )
                .repr()
        );
    }

    #[test]
    fn test_sub_1() {
        let a = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(66),
            utils::rev([0b1011_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let b = SmallFloat::from_parts(
            4,
            Sign::Neg,
            Exponent::Finite(2),
            utils::rev([0b1001_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let mut c = SmallFloat::<2>::zero(72);
        add_different_sign(&mut c, a.sign(), &a, &b, Round::ToNearest);
        assert!(
            c.repr()
                == SmallFloat::from_parts(
                    72,
                    Sign::Pos,
                    Exponent::Finite(66),
                    utils::rev([
                        0b1010_1111_1111_1111_1111_1111_1111_1111_1111_1111_1111_1111_1111_1111_1111_1111,
                        0b0111_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000
                    ]),
                )
                .repr()
        );
    }

    #[test]
    fn test_sub_2() {
        let a = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(68),
            utils::rev([0b1011_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let b = SmallFloat::from_parts(
            4,
            Sign::Neg,
            Exponent::Finite(2),
            utils::rev([0b1001_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let mut c = SmallFloat::<2>::zero(72);
        add_different_sign(&mut c, a.sign(), &a, &b, Round::ToNearest);
        assert!(
            c.repr()
                == SmallFloat::from_parts(
                    72,
                    Sign::Pos,
                    Exponent::Finite(68),
                    utils::rev([
                        0b1010_1111_1111_1111_1111_1111_1111_1111_1111_1111_1111_1111_1111_1111_1111_1111,
                        0b1101_1100_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000
                    ]),
                )
                .repr()
        );
    }

    #[test]
    fn test_sub_3() {
        let a = SmallFloat::from_parts(
            72,
            Sign::Pos,
            Exponent::Finite(2),
            utils::rev([
                0b1011_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000,
                0b1011_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000,
            ]),
        );
        let b = SmallFloat::from_parts(
            72,
            Sign::Neg,
            Exponent::Finite(2),
            utils::rev([
                0b1011_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000,
                0b1010_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000,
            ]),
        );
        let mut c = SmallFloat::<2>::zero(2);
        add_different_sign(&mut c, a.sign(), &a, &b, Round::ToNearest);
        assert!(
            c.repr()
                == SmallFloat::from_parts(
                    2,
                    Sign::Pos,
                    Exponent::Finite(2 - 67),
                    utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
                )
                .repr()
        );
    }
}
