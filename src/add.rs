use super::*;

pub fn add_same_sign(dst: &mut BigFloat, lhs_sign: Sign, lhs: &BigFloat, rhs: &BigFloat, rnd: Round) -> Approx {
    let (lhs, rhs) = if (lhs.sign_biased_exponent & consts::BIASED_EXPONENT_MASK) >= (rhs.sign_biased_exponent & consts::BIASED_EXPONENT_MASK) {
        (lhs, rhs)
    } else {
        (rhs, lhs)
    };

    let lhs_exp = lhs.exponent();
    let rhs_exp = rhs.exponent();
    if rhs_exp == Exponent::Zero || lhs_exp == Exponent::NaN || lhs_exp == Exponent::Inf {
        return dst.copy_from(lhs, rnd);
    }

    dst.sign_biased_exponent = (lhs.sign_biased_exponent & consts::BIASED_EXPONENT_MASK) | ((lhs_sign.is_negative() as u64) << 63);

    let rnd = rnd.with_sign(lhs.sign());
    let (Exponent::Finite(lhs_exp), Exponent::Finite(rhs_exp)) = (lhs_exp, rhs_exp) else {
        unreachable!()
    };

    let exp_shift = (lhs_exp as u64).wrapping_sub(rhs_exp as u64);
    let exp_large_shift = exp_shift / consts::LIMB_BITS;
    let exp_small_shift = exp_shift % consts::LIMB_BITS;

    let dst_prec = dst.precision_bits();

    let zero = consts::LIMB_ZERO;
    let one = consts::LIMB_ONE;

    let mut msb = false;
    let mut lsb_any_word = consts::LIMB_ZERO;
    let mut lsb_any = false;

    let dst_len = dst.mantissa_len() as u64;
    let dst_bits_mod_limb = dst_prec % consts::LIMB_BITS;

    let idx_that_contains_msb = if dst_bits_mod_limb == 0 { dst_len } else { dst_len - 1 };
    let msb_pos = consts::LIMB_BITS - dst_bits_mod_limb - 1;
    let ulp_pos = (consts::LIMB_BITS - dst_bits_mod_limb) % consts::LIMB_BITS;
    let ulp = one.shl(ulp_pos);

    let lhs_len = lhs.mantissa_len() as u64;
    let rhs_len = rhs.mantissa_len();
    let rhs_len_shifted = if exp_small_shift == 0 {
        rhs_len as u64 + exp_large_shift
    } else {
        rhs_len as u64 + exp_large_shift + 1
    };

    let max = Ord::max;
    let mut sum_iter = core::iter::from_fn({
        let mut i = max(dst_len, max(lhs_len, rhs_len_shifted));
        let mut carry = false;
        move || {
            if i == 0 {
                None
            } else {
                i -= 1;

                let l = if i < lhs_len {
                    let idx = lhs_len - i - 1;
                    lhs.mantissa()[idx as usize]
                } else {
                    zero
                };

                let r = if exp_small_shift == 0 {
                    if i < exp_large_shift || i >= rhs_len_shifted {
                        zero
                    } else {
                        let idx = rhs_len_shifted - i - 1;
                        rhs.mantissa()[idx as usize]
                    }
                } else {
                    let (r0, r1) = if i >= rhs_len_shifted {
                        (zero, zero)
                    } else if i + 1 == rhs_len_shifted {
                        (zero, rhs.mantissa()[0])
                    } else {
                        let idx0 = (rhs_len_shifted - i - 2) as usize;
                        let idx1 = (rhs_len_shifted - i - 1) as usize;

                        let r0 = if idx0 < rhs_len { rhs.mantissa()[idx0] } else { zero };
                        let r1 = if idx1 < rhs_len { rhs.mantissa()[idx1] } else { zero };
                        (r0, r1)
                    };

                    r0.shr(exp_small_shift) | r1.shl(consts::LIMB_BITS - exp_small_shift)
                };

                let (sum, carry_0) = l.overflowing_add(r);
                let (sum, carry_1) = sum.overflowing_add(carry as Limb);
                carry = carry_0 | carry_1;
                Some((sum, carry))
            }
        }
    });

    for i in (dst_len..max(dst_len, max(lhs_len, rhs_len_shifted))).rev() {
        let sum = sum_iter.next().unwrap().0;

        if i > idx_that_contains_msb {
            lsb_any_word |= sum;
        } else if i == idx_that_contains_msb {
            msb = (sum & one.shl(msb_pos)) == one.shl(msb_pos);
            lsb_any_word |= sum & one.shl(msb_pos).wrapping_sub(one);
            lsb_any = lsb_any_word != zero;
        }
    }

    let i = dst_len - 1;
    let (mut sum, mut carry) = sum_iter.next().unwrap();
    if i == idx_that_contains_msb {
        msb = (sum & one.shl(msb_pos)) == one.shl(msb_pos);
        lsb_any_word |= sum & one.shl(msb_pos).wrapping_sub(one);
        lsb_any = lsb_any_word != zero;
    }
    let mut result_lsb = (sum & ulp) == ulp;
    dst.mantissa_mut()[0] = sum;

    for i in (0..dst_len - 1).rev() {
        let idx = (dst_len - i - 1) as usize;
        (sum, carry) = sum_iter.next().unwrap();
        dst.mantissa_mut()[idx] = sum;
    }

    if carry {
        // shift to the right by 1
        lsb_any |= msb;
        msb = result_lsb;

        dst.sign_biased_exponent += 1;
        let dst_m = dst.mantissa_mut();
        for i in 0..dst_len as usize - 1 {
            dst_m[i] = dst_m[i].shr(1) | dst_m[i + 1].shl(consts::LIMB_BITS - 1);
        }
        let i = dst_len as usize - 1;
        dst_m[i] = dst_m[i].shr(1) | one.shl(consts::LIMB_BITS - 1);
        dst_m[0] &= !(one.shl(ulp_pos).wrapping_sub(one));
        result_lsb = (dst_m[0] & ulp) == ulp;
    }
    if dst.sign_biased_exponent & consts::BIASED_EXPONENT_MASK == consts::BIASED_EXPONENT_INF {
        return Approx::Overflow;
    }

    if rnd == RoundKnownSign::AwayFromZero || (rnd == RoundKnownSign::ToNearest && msb && (lsb_any || result_lsb)) {
        // add 1ulp
        let dst_m = dst.mantissa_mut();

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
            Approx::from_sign(lhs_sign)
        }
    } else {
        Approx::from_sign(lhs_sign.neg())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use equator::assert;

    #[test]
    fn test_add_same_sign_same_exp() {
        let a = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(4),
            utils::rev([0b1011_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let b = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(4),
            utils::rev([0b1011_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );

        let mut c = SmallFloat::from_parts(4, Sign::Neg, Exponent::Finite(2), utils::rev([u64::MAX]));
        assert!(add_same_sign(&mut c, a.sign(), &a, &b, Round::ToNearest) == Approx::Exact);

        assert!(
            c.repr()
                == SmallFloat::from_parts(
                    4,
                    Sign::Pos,
                    Exponent::Finite(5),
                    utils::rev([0b1011_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000])
                )
                .repr()
        );
    }

    #[test]
    fn test_add_same_sign_diff_exp() {
        //  10110
        //   1011
        // 100001
        let a = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(4),
            utils::rev([0b1011_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let b = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(5),
            utils::rev([0b1011_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );

        {
            let mut c = SmallFloat::from_parts(6, Sign::Neg, Exponent::Finite(2), utils::rev([u64::MAX]));
            assert!(add_same_sign(&mut c, a.sign(), &a, &b, Round::ToNearest) == Approx::Exact);
            assert!(
                c.repr()
                    == SmallFloat::from_parts(
                        6,
                        Sign::Pos,
                        Exponent::Finite(6),
                        utils::rev([0b1000_0100_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000])
                    )
                    .repr()
            );
        }
        {
            let mut c = SmallFloat::from_parts(5, Sign::Neg, Exponent::Finite(2), utils::rev([u64::MAX]));
            assert!(add_same_sign(&mut c, a.sign(), &a, &b, Round::ToNearest) == Approx::LessThanExact);
            assert!(
                c.repr()
                    == SmallFloat::from_parts(
                        5,
                        Sign::Pos,
                        Exponent::Finite(6),
                        utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000])
                    )
                    .repr()
            );
        }
        {
            let mut c = SmallFloat::from_parts(5, Sign::Neg, Exponent::Finite(2), utils::rev([u64::MAX]));
            assert!(add_same_sign(&mut c, a.sign(), &a, &b, Round::AwayFromZero) == Approx::GreaterThanExact);
            assert!(
                c.repr()
                    == SmallFloat::from_parts(
                        5,
                        Sign::Pos,
                        Exponent::Finite(6),
                        utils::rev([0b1000_1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000])
                    )
                    .repr()
            );
        }
        {
            let mut c = SmallFloat::from_parts(4, Sign::Neg, Exponent::Finite(2), utils::rev([u64::MAX]));
            assert!(add_same_sign(&mut c, a.sign(), &a, &b, Round::AwayFromZero) == Approx::GreaterThanExact);
            assert!(
                c.repr()
                    == SmallFloat::from_parts(
                        4,
                        Sign::Pos,
                        Exponent::Finite(6),
                        utils::rev([0b1001_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000])
                    )
                    .repr()
            );
        }
        {
            let mut c = SmallFloat::from_parts(4, Sign::Neg, Exponent::Finite(2), utils::rev([u64::MAX]));
            assert!(add_same_sign(&mut c, a.sign(), &a, &b, Round::ToNearest) == Approx::LessThanExact);
            assert!(
                c.repr()
                    == SmallFloat::from_parts(
                        4,
                        Sign::Pos,
                        Exponent::Finite(6),
                        utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000])
                    )
                    .repr()
            );
        }
    }

    #[test]
    fn test_add_same_sign_large_shift() {
        let a = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(68),
            utils::rev([0b1011_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let b = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(4),
            utils::rev([0b1011_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );

        let mut c = SmallFloat::from_parts(68, Sign::Neg, Exponent::Finite(2), utils::rev([u64::MAX; 2]));
        assert!(add_same_sign(&mut c, a.sign(), &a, &b, Round::ToNearest) == Approx::Exact);
        assert!(
            c.repr()
                == SmallFloat::from_parts(
                    68,
                    Sign::Pos,
                    Exponent::Finite(68),
                    utils::rev([
                        0b1011_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000,
                        0b1011_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000
                    ])
                )
                .repr()
        );

        let mut c = SmallFloat::from_parts(6, Sign::Neg, Exponent::Finite(2), utils::rev([u64::MAX; 2]));
        assert!(add_same_sign(&mut c, a.sign(), &a, &b, Round::ToNearest) == Approx::LessThanExact);
        assert!(
            c.repr()
                == SmallFloat::from_parts(
                    6,
                    Sign::Pos,
                    Exponent::Finite(68),
                    utils::rev([0b1011_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000,])
                )
                .repr()
        );

        let mut c = SmallFloat::from_parts(6, Sign::Neg, Exponent::Finite(2), utils::rev([u64::MAX; 2]));
        assert!(add_same_sign(&mut c, a.sign(), &a, &b, Round::Up) == Approx::GreaterThanExact);
        assert!(
            c.repr()
                == SmallFloat::from_parts(
                    6,
                    Sign::Pos,
                    Exponent::Finite(68),
                    utils::rev([0b1011_0100_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000,])
                )
                .repr()
        );

        let mut c = SmallFloat::from_parts(4, Sign::Neg, Exponent::Finite(2), utils::rev([u64::MAX; 2]));
        assert!(add_same_sign(&mut c, a.sign(), &a, &b, Round::Up) == Approx::GreaterThanExact);
        assert!(
            c.repr()
                == SmallFloat::from_parts(
                    4,
                    Sign::Pos,
                    Exponent::Finite(68),
                    utils::rev([0b1100_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000,])
                )
                .repr()
        );
    }

    #[test]
    fn test_add_same_sign_small_large_shift() {
        let a = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(68),
            utils::rev([0b1011_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let b = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(7),
            utils::rev([0b1011_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );

        let mut c = SmallFloat::from_parts(68, Sign::Neg, Exponent::Finite(2), utils::rev([u64::MAX; 2]));
        assert!(add_same_sign(&mut c, a.sign(), &a, &b, Round::ToNearest) == Approx::Exact);
        assert!(
            c.repr()
                == SmallFloat::from_parts(
                    68,
                    Sign::Pos,
                    Exponent::Finite(68),
                    utils::rev([
                        0b1011_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0101,
                        0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000
                    ])
                )
                .repr()
        );

        let mut c = SmallFloat::from_parts(6, Sign::Neg, Exponent::Finite(2), utils::rev([u64::MAX; 2]));
        assert!(add_same_sign(&mut c, a.sign(), &a, &b, Round::ToNearest) == Approx::LessThanExact);
        assert!(
            c.repr()
                == SmallFloat::from_parts(
                    6,
                    Sign::Pos,
                    Exponent::Finite(68),
                    utils::rev([0b1011_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000,])
                )
                .repr()
        );

        let mut c = SmallFloat::from_parts(6, Sign::Neg, Exponent::Finite(2), utils::rev([u64::MAX; 2]));
        assert!(add_same_sign(&mut c, a.sign(), &a, &b, Round::Up) == Approx::GreaterThanExact);
        assert!(
            c.repr()
                == SmallFloat::from_parts(
                    6,
                    Sign::Pos,
                    Exponent::Finite(68),
                    utils::rev([0b1011_0100_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000,])
                )
                .repr()
        );

        let mut c = SmallFloat::from_parts(4, Sign::Neg, Exponent::Finite(2), utils::rev([u64::MAX; 2]));
        assert!(add_same_sign(&mut c, a.sign(), &a, &b, Round::Up) == Approx::GreaterThanExact);
        assert!(
            c.repr()
                == SmallFloat::from_parts(
                    4,
                    Sign::Pos,
                    Exponent::Finite(68),
                    utils::rev([0b1100_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000,])
                )
                .repr()
        );
    }

    #[test]
    fn test_add_same_sign_small_very_large_shift() {
        let a = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(132),
            utils::rev([0b1011_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let b = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(7),
            utils::rev([0b1011_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );

        let mut c = SmallFloat::from_parts(132, Sign::Neg, Exponent::Finite(2), utils::rev([u64::MAX; 3]));
        assert!(add_same_sign(&mut c, a.sign(), &a, &b, Round::ToNearest) == Approx::Exact);
        assert!(
            c.repr()
                == SmallFloat::from_parts(
                    132,
                    Sign::Pos,
                    Exponent::Finite(132),
                    utils::rev([
                        0b1011_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000,
                        0b0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0101,
                        0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000,
                    ])
                )
                .repr()
        );
    }
}
