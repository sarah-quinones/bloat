use super::*;
use dyn_stack::PodStack;

mod schoolbook;

pub fn mul(dst: &mut BigFloat, lhs: &BigFloat, rhs: &BigFloat, rnd: Round, stack: PodStack<'_>) -> Approx {
    let sign = if lhs.sign().is_positive() { rhs.sign() } else { rhs.sign().neg() };
    let mut exp = match (lhs.exponent(), rhs.exponent()) {
        (Exponent::NaN, _) | (_, Exponent::NaN) | (Exponent::Zero, Exponent::Inf) | (Exponent::Inf, Exponent::Zero) => {
            dst.sign_biased_exponent = make_sign_and_biased_exponent(sign, Exponent::NaN);
            return Approx::Exact;
        }
        (Exponent::Finite(_), Exponent::Inf) | (Exponent::Inf, Exponent::Finite(_)) | (Exponent::Inf, Exponent::Inf) => {
            dst.sign_biased_exponent = make_sign_and_biased_exponent(sign, Exponent::Inf);
            return Approx::Exact;
        }
        (Exponent::Finite(_), Exponent::Zero) | (Exponent::Zero, Exponent::Finite(_)) | (Exponent::Zero, Exponent::Zero) => {
            dst.sign_biased_exponent = make_sign_and_biased_exponent(sign, Exponent::Zero);
            return Approx::Exact;
        }
        (Exponent::Finite(lhs_exp), Exponent::Finite(rhs_exp)) => lhs_exp + rhs_exp,
    };

    let lhs_m = lhs.mantissa();
    let rhs_m = rhs.mantissa();

    let (full_mul, _) = temp_big_float_uninit(consts::LIMB_BITS * (lhs_m.len() + rhs_m.len()) as u64, stack);

    schoolbook::mul_bigint(full_mul.full_mantissa_mut(), lhs_m, rhs_m);
    shl_1_if_needed(full_mul, &mut exp);

    if exp > consts::MAX_EXPONENT_INCLUSIVE {
        dst.sign_biased_exponent = make_sign_and_biased_exponent(sign, Exponent::Inf);
        return Approx::Exact;
    }
    if exp < consts::MIN_EXPONENT_INCLUSIVE {
        dst.sign_biased_exponent = make_sign_and_biased_exponent(sign, Exponent::Zero);
        return Approx::Exact;
    }

    full_mul.sign_biased_exponent = make_sign_and_biased_exponent(sign, Exponent::Finite(exp));
    dst.copy_from(&full_mul, rnd)
}

#[inline]
fn shl_1_if_needed(full_mul: &mut BigFloat, exp: &mut i64) {
    let full_mul = full_mul.full_mantissa_mut();
    let msb = consts::LIMB_ONE.shl(consts::LIMB_BITS - 1);
    if (full_mul[full_mul.len() - 1] & msb) != msb {
        for i in (1..full_mul.len()).rev() {
            full_mul[i] = full_mul[i].shl(1) | full_mul[i - 1].shr(consts::LIMB_BITS - 1)
        }
        full_mul[0] = full_mul[0].shl(1);
        *exp -= 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use equator::assert;

    #[test]
    fn test_mul_one_one() {
        let a = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(1),
            utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let b = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(1),
            utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );

        let mut c = SmallFloat::<1>::zero(4);
        assert!(mul(&mut c, &a, &b, Round::ToNearest, PodStack::new(bytemuck::cast_slice_mut(&mut [0u64; 4]))) == Approx::Exact);
        assert!(
            c.repr()
                == SmallFloat::from_parts(
                    4,
                    Sign::Pos,
                    Exponent::Finite(1),
                    utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
                )
                .repr()
        );
    }

    #[test]
    fn test_mul_one_half() {
        let a = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(1),
            utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let b = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(0),
            utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );

        let mut c = SmallFloat::<1>::zero(4);
        assert!(mul(&mut c, &a, &b, Round::ToNearest, PodStack::new(bytemuck::cast_slice_mut(&mut [0u64; 4]))) == Approx::Exact);
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
    fn test_mul_half_half() {
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
            utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );

        let mut c = SmallFloat::<1>::zero(4);
        assert!(mul(&mut c, &a, &b, Round::ToNearest, PodStack::new(bytemuck::cast_slice_mut(&mut [0u64; 4]))) == Approx::Exact);
        assert!(
            c.repr()
                == SmallFloat::from_parts(
                    4,
                    Sign::Pos,
                    Exponent::Finite(-1),
                    utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
                )
                .repr()
        );
    }

    #[test]
    fn test_mul_0_75() {
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
        assert!(mul(&mut c, &a, &b, Round::ToNearest, PodStack::new(bytemuck::cast_slice_mut(&mut [0u64; 4]))) == Approx::Exact);
        assert!(
            c.repr()
                == SmallFloat::from_parts(
                    4,
                    Sign::Pos,
                    Exponent::Finite(0),
                    utils::rev([0b1001_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
                )
                .repr()
        );
    }
}
