use super::*;

mod newton;

pub fn sqrt(dst: &mut BigFloat, x: &BigFloat, rnd: Round, stack: PodStack) -> Approx {
    let exp = match x.exponent() {
        Exponent::Zero => {
            dst.sign_biased_exponent = make_sign_and_biased_exponent(x.sign(), Exponent::Zero);
            return Approx::Exact;
        }
        Exponent::Finite(exp) => match x.sign() {
            Sign::Neg => {
                dst.sign_biased_exponent = make_sign_and_biased_exponent(Sign::Pos, Exponent::NaN);
                return Approx::Exact;
            }
            Sign::Pos => exp,
        },
        Exponent::Inf => {
            dst.sign_biased_exponent = match x.sign() {
                Sign::Neg => make_sign_and_biased_exponent(Sign::Pos, Exponent::NaN),
                Sign::Pos => make_sign_and_biased_exponent(Sign::Pos, Exponent::Inf),
            };
            return Approx::Exact;
        }
        Exponent::NaN => {
            dst.sign_biased_exponent = make_sign_and_biased_exponent(Sign::Pos, Exponent::NaN);
            return Approx::Exact;
        }
    };

    // we want 1 + max(dst.mantissa_len() + 1, x.mantissa_len() / 2) limbs

    let nlimbs = Ord::max(x.mantissa_len() / 2 + 1, dst.mantissa_len()) + 1;
    let (sqrt, mut stack) = temp_big_float_uninit(nlimbs as u64 * consts::LIMB_BITS, stack);
    sqrt.sign_biased_exponent = make_sign_and_biased_exponent(Sign::Pos, Exponent::Finite((exp + 1) / 2));
    let isqrt = sqrt.mantissa_mut();

    // we need nlimbs = (x.mantissa_len() * (LIMB_BITS / 2) + lshift.div_ceil(2)) / LIMB_BITS
    // nlimbs * LIMB_BITS - x.mantissa_len() * (LIMB_BITS / 2) = lshift.div_ceil(2)
    let lshift = (nlimbs as u64 * consts::LIMB_BITS - x.mantissa_len() as u64 * (consts::LIMB_BITS / 2)) * 2 - (exp as u64 % 2);
    newton::isqrt(isqrt, x.mantissa(), lshift, stack.rb_mut());

    let (isqrt2, mut stack) = stack.make_raw::<Limb>(nlimbs * 2);
    mul::isqr(isqrt2, isqrt, stack.rb_mut());

    // compare isqrt2 with x.mantissa() << lshift
    let large_lshift = (lshift / consts::LIMB_BITS) as usize;
    let small_lshift = lshift % consts::LIMB_BITS;
    let mut eq = true;
    for i in 0..large_lshift {
        if isqrt2[i] != consts::LIMB_ZERO {
            eq = false;
            break;
        }
    }

    if small_lshift == 0 {
        if eq {
            if isqrt2[large_lshift] != x.mantissa()[0].shl(small_lshift) {
                eq = false;
            }
        }

        if eq {
            for i in 1..x.mantissa_len() {
                if isqrt2[large_lshift + i] != x.mantissa()[i].shl(small_lshift) | x.mantissa()[i + 1].shr(consts::LIMB_BITS - small_lshift) {
                    eq = false;
                    break;
                }
            }
        }
    } else {
        if eq {
            for i in 0..x.mantissa_len() {
                if isqrt2[i + large_lshift] != x.mantissa()[i] {
                    eq = false;
                    break;
                }
            }
        }
    }

    if !eq {
        isqrt[0] |= consts::LIMB_ONE;
    }

    dst.copy_from(&sqrt, rnd)
}

#[cfg(test)]
mod tests {
    use super::*;
    use equator::assert;

    #[test]
    fn test_sqrt_0() {
        let x = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(2),
            utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );

        let mut y = SmallFloat::<2>::zero(99);
        assert!(sqrt(&mut y, &x, Round::Down, PodStack::new(bytemuck::cast_slice_mut(&mut [0u64; 100]))) == Approx::LessThanExact);
        assert!(
            y.repr()
                == SmallFloat::from_parts(
                    99,
                    Sign::Pos,
                    Exponent::Finite(1),
                    utils::rev([
                        0b1011010100000100111100110011001111111001110111100110010010000100,
                        0b0101100101111101100010011011001101100000000000000000000000000000,
                    ])
                )
                .repr()
        );
    }

    #[test]
    fn test_sqrt_1() {
        let x = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(3),
            utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );

        let mut y = SmallFloat::<2>::zero(99);
        assert!(sqrt(&mut y, &x, Round::Down, PodStack::new(bytemuck::cast_slice_mut(&mut [0u64; 100]))) == Approx::LessThanExact);
        assert!(
            y.repr()
                == SmallFloat::from_parts(
                    99,
                    Sign::Pos,
                    Exponent::Finite(2),
                    utils::rev([
                        0b1000000000000000000000000000000000000000000000000000000000000000,
                        0b0000000000000000000000000000000000000000000000000000000000000000,
                    ])
                )
                .repr()
        );
    }

    #[test]
    fn test_sqrt_2() {
        let x = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(2),
            utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );

        let mut mem = [0u64; 1000];
        let mut stack = PodStack::new(bytemuck::cast_slice_mut(&mut mem));

        let prec = 10000;

        // sqrt(2)^2 - 2
        let mut y = BigFloat::zero(prec);
        let mut y2 = BigFloat::zero(prec);
        let mut y2_minus_x = BigFloat::zero(prec);

        math::sqrt(&mut y, &x, Round::ToNearest, stack.rb_mut());
        math::mul(&mut y2, &y, &y, Round::ToNearest, stack.rb_mut());
        math::sub(&mut y2_minus_x, &y2, &x, Round::ToNearest);

        let as_f64 = y2_minus_x.to_f64(Round::ToNearest).0;
        let mantissa = y.mantissa();

        let prec = prec as u32;

        let x = rug::Float::with_val(prec, 2.0);
        let y = rug::Float::with_val(prec, x.sqrt_ref());
        let y2 = rug::Float::with_val(prec, &y * &y);
        let y2_minus_x = rug::Float::with_val(prec, &y2 - &x);

        let raw = unsafe { &*y.as_raw() };
        let mantissa_target = unsafe { core::slice::from_raw_parts(raw.d.as_ptr(), prec.div_ceil(u64::BITS) as usize) };

        assert!(as_f64 == y2_minus_x.to_f64());
        assert!(mantissa == mantissa_target);
    }
}
