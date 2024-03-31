use super::*;

mod schoolbook;

#[inline]
pub fn idivrem_normalized_req(lhs_len: usize, rhs_len: usize) -> Result<StackReq, SizeOverflow> {
    _ = lhs_len;
    _ = rhs_len;
    Ok(StackReq::empty())
}
#[inline]
pub fn idivrem_normalized(quo: &mut [Limb], rem: Option<&mut [Limb]>, lhs: &mut [Limb], rhs: &[Limb], stack: PodStack<'_>) {
    _ = stack;
    schoolbook::divrem_bigint(quo, rem, lhs, rhs)
}

pub fn div(dst: &mut BigFloat, lhs: &BigFloat, rhs: &BigFloat, rnd: Round, stack: PodStack<'_>) -> Approx {
    let sign = if lhs.sign() == rhs.sign() { Sign::Pos } else { Sign::Neg };
    let mut exp = match (lhs.exponent(), rhs.exponent()) {
        (Exponent::NaN, _) | (_, Exponent::NaN) | (Exponent::Zero, Exponent::Zero) | (Exponent::Inf, Exponent::Inf) => {
            dst.sign_biased_exponent = make_sign_and_biased_exponent(sign, Exponent::NaN);
            return Approx::Exact;
        }
        (Exponent::Finite(_), Exponent::Zero) | (Exponent::Inf, Exponent::Finite(_)) | (Exponent::Inf, Exponent::Zero) => {
            dst.sign_biased_exponent = make_sign_and_biased_exponent(sign, Exponent::Inf);
            return Approx::Exact;
        }
        (Exponent::Finite(_), Exponent::Inf) | (Exponent::Zero, Exponent::Finite(_)) | (Exponent::Zero, Exponent::Inf) => {
            dst.sign_biased_exponent = make_sign_and_biased_exponent(sign, Exponent::Zero);
            return Approx::Exact;
        }
        (Exponent::Finite(lhs_exp), Exponent::Finite(rhs_exp)) => lhs_exp - rhs_exp,
    };

    let lhs_m = lhs.mantissa();
    let rhs_m = rhs.mantissa();

    let quo_len = Ord::max(dst.mantissa_len() + 1, lhs_m.len()) + 1;

    let (lhs_bigint, stack) = stack.make_raw::<Limb>(quo_len + rhs_m.len());

    let len = lhs_bigint.len();
    lhs_bigint[..len - lhs_m.len() - 1].fill(consts::LIMB_ZERO);
    lhs_bigint[len - lhs_m.len() - 1..len - 1].copy_from_slice(lhs_m);
    lhs_bigint[len - 1] = consts::LIMB_ZERO;

    let (quo, stack) = temp_big_float_uninit(quo_len as u64 * consts::LIMB_BITS, stack);
    let quo_bigint = quo.mantissa_mut();
    let (rem_bigint, _) = stack.make_raw::<Limb>(rhs_m.len());

    schoolbook::divrem_bigint(quo_bigint, Some(rem_bigint), lhs_bigint, rhs_m);
    let need_shift = quo_bigint[quo_bigint.len() - 1] == 1;

    if need_shift {
        let quo = quo.mantissa_mut();
        for i in (1..quo.len()).rev() {
            quo[i] = quo[i].shl(consts::LIMB_BITS - 1) | quo[i - 1].shr(1);
        }
        quo[0] = quo[0].shl(consts::LIMB_BITS - 1);
        exp += 1;
    } else {
        let quo = quo.mantissa_mut();
        quo.copy_within(0..quo_len - 1, 1);
        quo[0] = consts::LIMB_ZERO;
    }

    if exp > consts::MAX_EXPONENT_INCLUSIVE {
        dst.sign_biased_exponent = make_sign_and_biased_exponent(sign, Exponent::Inf);
        return Approx::Overflow;
    }
    if exp < consts::MIN_EXPONENT_INCLUSIVE {
        dst.sign_biased_exponent = make_sign_and_biased_exponent(sign, Exponent::Zero);
        return Approx::Underflow;
    }

    dst.sign_biased_exponent = make_sign_and_biased_exponent(sign, Exponent::Finite(exp));
    quo.sign_biased_exponent = make_sign_and_biased_exponent(sign, Exponent::Finite(exp));

    let mut zero_rem = true;
    for &x in &*rem_bigint {
        if x != consts::LIMB_ZERO {
            zero_rem = false;
        }
    }

    if !zero_rem {
        // add a fake bit to influence the rounding the way we want
        quo.mantissa_mut()[0] |= consts::LIMB_ONE;
    }
    dst.copy_from(quo, rnd)
}

#[cfg(test)]
mod tests {
    use super::*;
    use equator::assert;

    #[test]
    fn test_div_0() {
        let a = SmallFloat::from_parts(4, Sign::Pos, Exponent::Finite(2), [consts::LIMB_HIGH_BIT]);
        let b = SmallFloat::from_parts(4, Sign::Pos, Exponent::Finite(2), [consts::LIMB_HIGH_BIT]);

        let mut c = SmallFloat::<1>::zero(4);
        div(
            &mut c,
            &a,
            &b,
            Round::ToNearest,
            PodStack::new(bytemuck::cast_slice_mut(&mut [0u64; 100])),
        );
        assert!(c.repr() == SmallFloat::from_parts(4, Sign::Pos, Exponent::Finite(1), [consts::LIMB_HIGH_BIT]).repr());
    }

    #[test]
    fn test_div_1() {
        let a = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(0),
            utils::rev([0b1100_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let b = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(0),
            utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );

        let mut c = SmallFloat::<1>::zero(4);
        div(
            &mut c,
            &a,
            &b,
            Round::ToNearest,
            PodStack::new(bytemuck::cast_slice_mut(&mut [0u64; 100])),
        );
        assert!(
            c.repr()
                == SmallFloat::from_parts(
                    4,
                    Sign::Pos,
                    Exponent::Finite(1),
                    utils::rev([0b1100_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000])
                )
                .repr()
        );
    }

    #[test]
    fn test_div_2() {
        let a = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(0),
            utils::rev([0b1100_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let b = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(0),
            utils::rev([0b1100_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );

        let mut c = SmallFloat::<1>::zero(4);
        div(
            &mut c,
            &a,
            &b,
            Round::ToNearest,
            PodStack::new(bytemuck::cast_slice_mut(&mut [0u64; 100])),
        );
        assert!(
            c.repr()
                == SmallFloat::from_parts(
                    4,
                    Sign::Pos,
                    Exponent::Finite(1),
                    utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000])
                )
                .repr()
        );
    }

    #[test]
    fn test_div_3() {
        let a = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(0),
            utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let b = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(0),
            utils::rev([0b1100_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );

        let mut c = SmallFloat::<1>::zero(16);
        assert!(
            div(
                &mut c,
                &a,
                &b,
                Round::ToNearest,
                PodStack::new(bytemuck::cast_slice_mut(&mut [0u64; 100])),
            ) == Approx::GreaterThanExact
        );
        // 0.666..
        assert!(
            c.repr()
                == SmallFloat::from_parts(
                    16,
                    Sign::Pos,
                    Exponent::Finite(0),
                    utils::rev([0b1010_1010_1010_1011_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000])
                )
                .repr()
        );
    }

    #[test]
    fn test_div_4() {
        let a = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(0),
            utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let b = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(0),
            utils::rev([0b1100_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );

        let mut c = SmallFloat::<1>::zero(15);
        assert!(
            div(
                &mut c,
                &a,
                &b,
                Round::ToNearest,
                PodStack::new(bytemuck::cast_slice_mut(&mut [0u64; 100])),
            ) == Approx::LessThanExact
        );
        // 0.666..
        assert!(
            c.repr()
                == SmallFloat::from_parts(
                    15,
                    Sign::Pos,
                    Exponent::Finite(0),
                    utils::rev([0b1010_1010_1010_1010_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000])
                )
                .repr()
        );
    }

    #[test]
    fn test_div_5() {
        let a = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(0),
            utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let b = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(0),
            utils::rev([0b1100_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );

        let mut c = SmallFloat::<1>::zero(15);
        assert!(
            div(
                &mut c,
                &a,
                &b,
                Round::AwayFromZero,
                PodStack::new(bytemuck::cast_slice_mut(&mut [0u64; 100])),
            ) == Approx::GreaterThanExact
        );
        // 0.666..
        assert!(
            c.repr()
                == SmallFloat::from_parts(
                    15,
                    Sign::Pos,
                    Exponent::Finite(0),
                    utils::rev([0b1010_1010_1010_1100_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000])
                )
                .repr()
        );
    }

    #[test]
    fn test_div_6() {
        let a = SmallFloat::from_parts(
            132,
            Sign::Pos,
            Exponent::Finite(0),
            utils::rev([
                0b1000_0001_0001_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000,
                0b0001_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000,
                0b0001_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000,
            ]),
        );
        let b = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(0),
            utils::rev([0b1100_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );

        let mut c = SmallFloat::<1>::zero(4);
        assert!(
            div(
                &mut c,
                &a,
                &b,
                Round::ToNearest,
                PodStack::new(bytemuck::cast_slice_mut(&mut [0u64; 100])),
            ) == Approx::GreaterThanExact
        );
        // 0.666..
        assert!(
            c.repr()
                == SmallFloat::from_parts(
                    4,
                    Sign::Pos,
                    Exponent::Finite(0),
                    utils::rev([0b1011_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000])
                )
                .repr()
        );
    }

    #[test]
    fn test_div_7() {
        let a = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(0),
            utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let b = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(0),
            utils::rev([0b1100_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );

        let mut c = SmallFloat::<2>::zero(68);
        assert!(
            div(
                &mut c,
                &a,
                &b,
                Round::ToNearest,
                PodStack::new(bytemuck::cast_slice_mut(&mut [0u64; 100])),
            ) == Approx::GreaterThanExact
        );
        // 0.666..
        assert!(
            c.repr()
                == SmallFloat::from_parts(
                    68,
                    Sign::Pos,
                    Exponent::Finite(0),
                    utils::rev([
                        0b1010_1010_1010_1010_1010_1010_1010_1010_1010_1010_1010_1010_1010_1010_1010_1010,
                        0b1011_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000
                    ])
                )
                .repr()
        );
    }
}
