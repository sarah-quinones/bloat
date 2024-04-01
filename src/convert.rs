use super::*;

pub fn to_f64(x: &BigFloat, rnd: Round) -> (f64, Approx) {
    let mut as_f64 = SmallFloat::<1>::zero(53);
    let approx = as_f64.copy_from(x, rnd);
    let mut f = (as_f64.sign().is_negative() as u64) << 63;
    match as_f64.exponent() {
        Exponent::Zero => {}
        Exponent::Finite(exp) => {
            let exp = exp - 1;
            if exp < -1022 {
                // zero
            } else if exp > 1023 {
                f |= 2047 << 52;
            } else {
                let exp = (exp + 1023) as u64;
                f |= exp << 52;
                f |= (as_f64.mantissa()[0] >> (64 - 53)) & ((1 << 52) - 1);
            }
        }
        Exponent::Inf => {
            f |= 2047 << 52;
        }
        Exponent::NaN => {
            f |= 2047 << 52;
            f |= (1 << 52) - 1;
        }
    }

    (f64::from_bits(f), approx)
}

pub fn from_f64(dst: &mut BigFloat, x: f64, rnd: Round) -> Approx {
    let sign = if x.is_sign_negative() { Sign::Neg } else { Sign::Pos };

    if x.is_nan() {
        dst.sign_biased_exponent = make_sign_and_biased_exponent(sign, Exponent::NaN);
        Approx::Exact
    } else if x.is_infinite() {
        dst.sign_biased_exponent = make_sign_and_biased_exponent(sign, Exponent::Inf);
        Approx::Exact
    } else if x == 0.0 {
        dst.sign_biased_exponent = make_sign_and_biased_exponent(sign, Exponent::Zero);
        Approx::Exact
    } else {
        let mut tmp = SmallFloat::<1>::zero(53);
        let x = x.to_bits();
        let biased_exp = (x & ((1u64 << 63) - 1)) >> 52;
        let exp = biased_exp as i64 - 1023;
        let mantissa = x & ((1u64 << 52) - 1);

        if biased_exp == 0 {
            let mantissa = mantissa << 11;
            let exp = exp - mantissa.leading_zeros() as i64 + 1;
            let mantissa = mantissa << mantissa.leading_zeros();
            tmp.sign_biased_exponent = make_sign_and_biased_exponent(sign, Exponent::Finite(exp + 1));
            tmp.mantissa_mut()[0] = mantissa;
        } else {
            let mantissa = ((1u64 << 52) | mantissa) << 11;
            tmp.sign_biased_exponent = make_sign_and_biased_exponent(sign, Exponent::Finite(exp + 1));
            tmp.mantissa_mut()[0] = mantissa;
        }

        dst.copy_from(&tmp, rnd)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use equator::assert;

    #[test]
    fn test_normal() {
        let x = 0.3;
        let mut y = SmallFloat::<1>::zero(53);
        assert!(from_f64(&mut y, x, Round::ToNearest) == Approx::Exact);
        assert!(to_f64(&y, Round::ToNearest).0 == x);
    }

    #[test]
    fn test_subnormal() {
        let x = f64::MIN_POSITIVE;
        let mut y = SmallFloat::<1>::zero(53);
        let mut z = SmallFloat::<1>::zero(53);
        let mut w = SmallFloat::<1>::zero(53);

        assert!(from_f64(&mut y, x, Round::ToNearest) == Approx::Exact);
        assert!(from_f64(&mut z, x / 2.0, Round::ToNearest) == Approx::Exact);
        assert!(from_f64(&mut w, x / 4.0, Round::ToNearest) == Approx::Exact);

        let Exponent::Finite(z_exp) = z.exponent() else { panic!() };
        z.sign_biased_exponent = make_sign_and_biased_exponent(z.sign(), Exponent::Finite(z_exp + 1));
        let Exponent::Finite(w_exp) = w.exponent() else { panic!() };
        w.sign_biased_exponent = make_sign_and_biased_exponent(w.sign(), Exponent::Finite(w_exp + 2));

        assert!(z.repr() == y.repr());
        assert!(w.repr() == y.repr());
    }
}
