use super::*;

mod schoolbook;

#[inline]
pub fn imul_req(lhs_len: usize, rhs_len: usize) -> Result<StackReq, SizeOverflow> {
    _ = lhs_len;
    _ = rhs_len;
    Ok(StackReq::empty())
}

#[inline]
pub fn isqr_req(input_len: usize) -> Result<StackReq, SizeOverflow> {
    _ = input_len;
    Ok(StackReq::empty())
}

#[inline]
pub fn imul(full_mul: &mut [Limb], lhs: &[Limb], rhs: &[Limb], stack: PodStack<'_>) {
    _ = stack;
    schoolbook::mul_bigint(full_mul, lhs, rhs);
}

#[inline]
pub fn isqr(full_mul: &mut [Limb], input: &[Limb], stack: PodStack<'_>) {
    _ = stack;
    schoolbook::mul_bigint(full_mul, input, input);
}

enum NegateKind {
    Negate,
    NoNegate,
}
enum MulKind {
    MulAdd,
    MulSub,
}

fn mul_generic(
    dst: &mut BigFloat,
    lhs: &BigFloat,
    rhs: &BigFloat,
    acc: Option<&BigFloat>,
    mul: MulKind,
    negate: NegateKind,
    rnd: Round,
    stack: PodStack<'_>,
) -> Approx {
    let sign = if lhs.sign() == rhs.sign() { Sign::Pos } else { Sign::Neg };
    let sign = match negate {
        NegateKind::Negate => sign.neg(),
        NegateKind::NoNegate => sign,
    };

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
            let mut zero = SmallFloat::<1>::zero(2);
            zero.sign_biased_exponent = make_sign_and_biased_exponent(sign, Exponent::Zero);

            return if let Some(acc) = acc {
                match mul {
                    MulKind::MulAdd => math::add(dst, &zero, acc, rnd),
                    MulKind::MulSub => math::sub(dst, &zero, acc, rnd),
                }
            } else {
                dst.sign_biased_exponent = make_sign_and_biased_exponent(sign, Exponent::Zero);
                Approx::Exact
            };
        }
        (Exponent::Finite(lhs_exp), Exponent::Finite(rhs_exp)) => lhs_exp + rhs_exp,
    };

    let lhs_m = lhs.mantissa();
    let rhs_m = rhs.mantissa();

    let (full_mul, _) = temp_big_float_uninit(consts::LIMB_BITS * (lhs_m.len() + rhs_m.len()) as u64, stack);

    schoolbook::mul_bigint(full_mul.full_mantissa_mut(), lhs_m, rhs_m);
    shl1_if_needed(full_mul, &mut exp);

    if exp > consts::MAX_EXPONENT_INCLUSIVE {
        dst.sign_biased_exponent = make_sign_and_biased_exponent(sign, Exponent::Inf);
        return Approx::Overflow;
    }
    if exp < consts::MIN_EXPONENT_INCLUSIVE {
        dst.sign_biased_exponent = make_sign_and_biased_exponent(sign, Exponent::Zero);
        return Approx::Underflow;
    }

    full_mul.sign_biased_exponent = make_sign_and_biased_exponent(sign, Exponent::Finite(exp));
    if let Some(acc) = acc {
        match mul {
            MulKind::MulAdd => math::add(dst, full_mul, acc, rnd),
            MulKind::MulSub => math::sub(dst, full_mul, acc, rnd),
        }
    } else {
        dst.copy_from(&full_mul, rnd)
    }
}

pub fn mul(dst: &mut BigFloat, lhs: &BigFloat, rhs: &BigFloat, rnd: Round, stack: PodStack<'_>) -> Approx {
    mul_generic(dst, lhs, rhs, None, MulKind::MulAdd, NegateKind::NoNegate, rnd, stack)
}

pub fn mul_add(dst: &mut BigFloat, lhs: &BigFloat, rhs: &BigFloat, acc: &BigFloat, rnd: Round, stack: PodStack<'_>) -> Approx {
    mul_generic(dst, lhs, rhs, Some(acc), MulKind::MulAdd, NegateKind::NoNegate, rnd, stack)
}

pub fn mul_sub(dst: &mut BigFloat, lhs: &BigFloat, rhs: &BigFloat, acc: &BigFloat, rnd: Round, stack: PodStack<'_>) -> Approx {
    mul_generic(dst, lhs, rhs, Some(acc), MulKind::MulSub, NegateKind::NoNegate, rnd, stack)
}

pub fn negate_mul_add(dst: &mut BigFloat, lhs: &BigFloat, rhs: &BigFloat, acc: &BigFloat, rnd: Round, stack: PodStack<'_>) -> Approx {
    mul_generic(dst, lhs, rhs, Some(acc), MulKind::MulAdd, NegateKind::Negate, rnd, stack)
}

pub fn negate_mul_sub(dst: &mut BigFloat, lhs: &BigFloat, rhs: &BigFloat, acc: &BigFloat, rnd: Round, stack: PodStack<'_>) -> Approx {
    mul_generic(dst, lhs, rhs, Some(acc), MulKind::MulSub, NegateKind::Negate, rnd, stack)
}

#[inline]
fn shl1_if_needed(full_mul: &mut BigFloat, exp: &mut i64) {
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
