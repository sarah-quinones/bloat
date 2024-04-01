use super::*;

fn remquo_finite_direct_req(rem_prec: u64, quo_len: usize, lhs_prec: u64, rhs_prec: u64) -> Result<StackReq, SizeOverflow> {
    let fquo_prec = (quo_len as u64).checked_mul(consts::LIMB_BITS).ok_or(SizeOverflow)?;

    StackReq::try_all_of([
        temp_big_float_req(fquo_prec)?,
        StackReq::try_any_of([
            math::div_req(fquo_prec, lhs_prec, rhs_prec)?,
            math::mul_req(rem_prec, fquo_prec, rhs_prec)?,
        ])?,
    ])
}
fn remquo_finite_req(rem_prec: u64, quo_len: usize, lhs_prec: u64, rhs_prec: u64) -> Result<StackReq, SizeOverflow> {
    let y_len = (rhs_prec.div_ceil(consts::LIMB_BITS)) as usize;

    let wide = StackReq::try_new::<Limb>(y_len * 2 + 1)?;
    let pow2_mod_y = StackReq::try_new::<Limb>(y_len)?;
    let __quo = StackReq::try_new::<Limb>(y_len)?;
    let x_pow2_mod_y = StackReq::try_new::<Limb>(y_len)?;

    let idiv = div::idivrem_normalized_req(y_len * 2 + 1, y_len)?;
    let imul = mul::imul_req(y_len, y_len)?;

    let frem = temp_big_float_req((y_len as u64 + 1).checked_mul(consts::LIMB_BITS).ok_or(SizeOverflow)?)?;

    StackReq::try_all_of([
        pow2_mod_y,
        StackReq::try_any_of([
            StackReq::try_all_of([frem, remquo_finite_direct_req(rem_prec, quo_len, lhs_prec, rhs_prec)?])?,
            StackReq::try_all_of([wide, __quo, x_pow2_mod_y, StackReq::try_any_of([idiv, imul])?])?,
        ])?,
    ])
}

pub fn remquo_req(rem_prec: u64, quo_len: usize, lhs_prec: u64, rhs_prec: u64) -> Result<StackReq, SizeOverflow> {
    remquo_finite_req(rem_prec, quo_len, lhs_prec, rhs_prec)
}

fn remquo_finite_direct(rem: &mut BigFloat, quo: &mut [Limb], lhs: &BigFloat, rhs: &BigFloat, rnd: Round, stack: PodStack<'_>) -> Approx {
    let Exponent::Finite(lhs_exp) = lhs.exponent() else { panic!() };
    let Exponent::Finite(rhs_exp) = rhs.exponent() else { panic!() };

    if lhs_exp <= rhs_exp - 1 {
        match (lhs.sign(), rhs.sign()) {
            (Sign::Neg, Sign::Neg) => {
                quo.fill(consts::LIMB_ONE);
                math::sub(rem, lhs, rhs, rnd)
            }
            (Sign::Neg, Sign::Pos) => {
                quo.fill(consts::LIMB_ONE);
                math::add(rem, lhs, rhs, rnd)
            }
            (Sign::Pos, _) => {
                quo.fill(consts::LIMB_ZERO);
                math::copy(rem, lhs, rnd)
            }
        }
    } else {
        let (fquo, mut stack) = temp_big_float_uninit((lhs_exp - rhs_exp + 2) as u64, stack);
        math::div(fquo, lhs, rhs, Round::Down, stack.rb_mut());
        fquo.sign_biased_exponent = make_sign_and_biased_exponent(Sign::Pos, fquo.exponent());

        let Exponent::Finite(fquo_exp) = fquo.exponent() else { panic!() };

        if fquo_exp == -1 {
            match (lhs.sign(), rhs.sign()) {
                (Sign::Neg, Sign::Neg) => {
                    quo.fill(consts::LIMB_ONE);
                    math::sub(rem, lhs, rhs, rnd)
                }
                (Sign::Neg, Sign::Pos) => {
                    quo.fill(consts::LIMB_ONE);
                    math::add(rem, lhs, rhs, rnd)
                }
                (Sign::Pos, _) => {
                    quo.fill(consts::LIMB_ZERO);
                    math::copy(rem, lhs, rnd)
                }
            }
        } else {
            // zero the least significant bit
            let fquo_bits_mod_limb = fquo.precision_bits() % consts::LIMB_BITS;
            let ulp_pos = (consts::LIMB_BITS - fquo_bits_mod_limb) % consts::LIMB_BITS;
            let ulp = consts::LIMB_ONE.shl(ulp_pos);
            let mut begin = 0;
            fquo.mantissa_mut()[begin] &= !ulp;
            if ulp_pos == consts::LIMB_BITS - 1 {
                begin += 1;
            }

            let ulp_pos = if fquo_exp == lhs_exp - rhs_exp {
                // zero the second least significant bit
                let ulp_pos = (ulp_pos + 1) % consts::LIMB_BITS;
                let ulp = consts::LIMB_ONE.shl(ulp_pos);

                fquo.mantissa_mut()[begin] &= !ulp;

                if ulp_pos == consts::LIMB_BITS - 1 {
                    begin += 1;
                }

                ulp_pos
            } else {
                ulp_pos
            };
            let ulp_pos = (ulp_pos + 1) % consts::LIMB_BITS;
            let ulp = consts::LIMB_ONE.shl(ulp_pos);

            {
                // copy fquo into quo
                let quo_len = quo.len();
                let fquo = &fquo.mantissa()[begin..];

                let zero = &consts::LIMB_ZERO;
                if ulp_pos == 0 {
                    for i in 0..quo_len {
                        quo[i] = *fquo.get(begin + i).unwrap_or(zero);
                    }
                } else {
                    for i in 0..quo_len {
                        let q0 = fquo.get(i).unwrap_or(zero).shr(ulp_pos);
                        let q1 = fquo.get(i + 1).unwrap_or(zero).shl(consts::LIMB_BITS - ulp_pos);
                        quo[i] = q0 | q1;
                    }
                }
            }

            if lhs.sign().is_negative() {
                for q in quo {
                    *q = !*q;
                }

                // fquo += 1;
                let mut carry;
                {
                    let fquo = &mut fquo.mantissa_mut()[begin..];

                    (fquo[0], carry) = fquo[0].overflowing_add(ulp);
                    for i in 1..fquo.len() {
                        (fquo[i], carry) = fquo[i].overflowing_add(carry as Limb);
                    }
                }
                if carry {
                    *fquo.mantissa_mut().last_mut().unwrap() = consts::LIMB_HIGH_BIT;
                    fquo.sign_biased_exponent = make_sign_and_biased_exponent(fquo.sign(), Exponent::Finite(fquo_exp + 1));
                }
            }

            if lhs.sign() == rhs.sign() {
                math::negate_mul_add(rem, fquo, rhs, lhs, rnd, stack)
            } else {
                math::mul_add(rem, fquo, rhs, lhs, rnd, stack)
            }
        }
    }
}

fn remquo_finite(rem: &mut BigFloat, quo: &mut [Limb], lhs: &BigFloat, rhs: &BigFloat, rnd: Round, stack: PodStack<'_>) -> Approx {
    let x = lhs;
    let y = rhs;

    let Exponent::Finite(x_exp) = x.exponent() else { panic!() };
    let Exponent::Finite(y_exp) = y.exponent() else { panic!() };
    let x_bitlen = x.mantissa_len() as u64 * consts::LIMB_BITS;
    let y_bitlen = y.mantissa_len() as u64 * consts::LIMB_BITS;
    let q_bitlen = quo.len() as u64 * consts::LIMB_BITS;

    // y = 0.ym * 2^y_exp
    // y * 2^(BITS * y.mantissa_len) = ym * 2^y_exp

    // y * 2^py
    // y * 2^(BITS * y.mantissa_len() - y_exp)
    let ym = y.mantissa();

    // x * 2^px
    // x * 2^(BITS * x.mantissa_len() - x_exp)
    // x * 2^py * 2^(BITS * (x.mantissa_len() - y.mantissa_len()) - (x_exp - y_exp))
    // i.e., x * 2^py = xm * 2^(BITS * (y.mantissa_len() - x.mantissa_len()) - (y_exp - x_exp))
    let xm = x.mantissa();

    if (y_bitlen.wrapping_sub(x_bitlen) as i64) >= y_exp - x_exp + q_bitlen as i64 {
        // xl        = 0
        // xh * 2^py = x * 2^p
        let p = (y_bitlen.wrapping_sub(x_bitlen) as i64).wrapping_sub(y_exp.wrapping_sub(x_exp).wrapping_add(q_bitlen as i64)) as u64;

        // https://en.wikipedia.org/wiki/Exponentiation_by_squaring
        let (y_pow2_mod_y, mut stack) = stack.make_raw::<Limb>(ym.len());
        let mut y_bit_pos = Some(0u64);
        {
            let (wide, stack) = stack.rb_mut().make_raw::<Limb>(ym.len() * 2 + 1);
            let (__quo, mut stack) = stack.make_raw::<Limb>(ym.len());
            let (x_pow2_mod_y, mut stack) = stack.rb_mut().make_raw::<Limb>(ym.len());

            wide.fill(consts::LIMB_ZERO);
            y_pow2_mod_y.fill(consts::LIMB_ZERO);
            x_pow2_mod_y.fill(consts::LIMB_ZERO);

            let mut x_bit_pos = Some(1u64);

            for idx in 0..u64::BITS as u64 {
                if p >> idx == 0 {
                    break;
                }
                let bit = ((p >> idx) & 1) == 1;
                if bit {
                    // y = x * y
                    if let Some(old) = y_bit_pos {
                        if let Some(x_bit_pos) = x_bit_pos {
                            let new = x_bit_pos + old;
                            if new < y_bitlen {
                                y_bit_pos = Some(new);
                            } else {
                                y_bit_pos = None;
                                wide.fill(consts::LIMB_ZERO);
                                wide[(new / consts::LIMB_BITS) as usize] = consts::LIMB_ONE.shl(new % consts::LIMB_BITS);
                                div::idivrem_normalized(__quo, Some(y_pow2_mod_y), wide, ym, stack.rb_mut());
                            }
                        } else {
                            y_pow2_mod_y[(old / consts::LIMB_BITS) as usize] = consts::LIMB_ONE.shl(old % consts::LIMB_BITS);
                            wide[ym.len() * 2] = consts::LIMB_ZERO;
                            mul::imul(&mut wide[..ym.len() * 2], &*x_pow2_mod_y, &*y_pow2_mod_y, stack.rb_mut());
                            div::idivrem_normalized(__quo, Some(y_pow2_mod_y), wide, ym, stack.rb_mut());
                        }
                    } else {
                        wide[ym.len() * 2] = consts::LIMB_ZERO;
                        mul::imul(&mut wide[..ym.len() * 2], &*x_pow2_mod_y, &*y_pow2_mod_y, stack.rb_mut());
                        div::idivrem_normalized(__quo, Some(y_pow2_mod_y), wide, ym, stack.rb_mut());
                    }
                }

                // x = x * x
                if let Some(old) = x_bit_pos {
                    let new = old * 2;
                    if new < y_bitlen {
                        x_bit_pos = Some(new);
                    } else {
                        x_bit_pos = None;
                        wide.fill(consts::LIMB_ZERO);
                        wide[(new / consts::LIMB_BITS) as usize] = consts::LIMB_ONE.shl(new % consts::LIMB_BITS);
                        div::idivrem_normalized(__quo, Some(x_pow2_mod_y), wide, ym, stack.rb_mut());
                    }
                } else {
                    wide[ym.len() * 2] = consts::LIMB_ZERO;
                    mul::isqr(&mut wide[..ym.len() * 2], &*x_pow2_mod_y, stack.rb_mut());
                    div::idivrem_normalized(__quo, Some(x_pow2_mod_y), wide, ym, stack.rb_mut());
                }
            }
        }

        let pow2_mod_y = y_pow2_mod_y;

        if let Some(y_bit_pos) = y_bit_pos {
            let (wide, stack) = stack.rb_mut().make_raw::<Limb>(ym.len() * 2 + 1);
            let (__quo, mut stack) = stack.make_raw::<Limb>(ym.len());
            pow2_mod_y[(y_bit_pos / consts::LIMB_BITS) as usize] = consts::LIMB_ONE.shl(y_bit_pos % consts::LIMB_BITS);
            wide[ym.len() * 2] = consts::LIMB_ZERO;
            mul::imul(&mut wide[..ym.len() * 2], &*xm, &*pow2_mod_y, stack.rb_mut());
            div::idivrem_normalized(__quo, Some(pow2_mod_y), wide, ym, stack.rb_mut());
        }

        let r = &*pow2_mod_y;

        // convert r from a bigint into a bigfloat
        let mut n = r.len();
        while n > 0 {
            n -= 1;

            if r[n] != consts::LIMB_ZERO {
                break;
            }
        }

        let lshift = r[n].leading_zeros() as u64;
        if lshift == consts::LIMB_BITS {
            // r is zero
            quo.fill(consts::LIMB_ZERO);
            rem.mantissa_mut().fill(consts::LIMB_ZERO);
            rem.sign_biased_exponent = make_sign_and_biased_exponent(Sign::Pos, Exponent::Zero);
            Approx::Exact
        } else {
            let (frem, stack) = temp_big_float_uninit((n + 1) as u64 * consts::LIMB_BITS, stack);
            let offset = (y_bitlen as i64).wrapping_sub(y_exp).wrapping_sub(q_bitlen as i64);
            frem.sign_biased_exponent = make_sign_and_biased_exponent(
                lhs.sign(),
                Exponent::Finite((n as u64 * consts::LIMB_BITS + (consts::LIMB_BITS - lshift)) as i64 - offset),
            );
            let frem_m = frem.mantissa_mut();
            frem_m[0] = r[0].shl(lshift);

            if lshift == 0 {
                for i in 1..n {
                    frem_m[i] = r[i].shl(lshift);
                }
            } else {
                for i in 1..n {
                    frem_m[i] = r[i].shl(lshift) | r[i - 1].shr(consts::LIMB_BITS - lshift);
                }
            }

            remquo_finite_direct(rem, quo, frem, rhs, rnd, stack)
        }
    } else {
        // (y_bitlen.wrapping_sub(x_bitlen) as i64) < y_exp - x_exp + q_bitlen as i64
        // x_exp - y_exp < q_bitlen + x_bitlen - y_bitlen

        remquo_finite_direct(rem, quo, lhs, rhs, rnd, stack)
    }
}

pub fn remquo(rem: &mut BigFloat, quo: &mut [Limb], lhs: &BigFloat, rhs: &BigFloat, rnd: Round, stack: PodStack<'_>) -> Approx {
    match (lhs.exponent(), rhs.exponent()) {
        (Exponent::NaN | Exponent::Inf, _) | (_, Exponent::NaN | Exponent::Zero) => {
            quo.fill(consts::LIMB_ZERO);
            rem.sign_biased_exponent = make_sign_and_biased_exponent(Sign::Pos, Exponent::NaN);
            Approx::Exact
        }
        (Exponent::Zero, _) | (Exponent::Finite(_), Exponent::Inf) => {
            quo.fill(consts::LIMB_ZERO);
            rem.copy_from(lhs, rnd)
        }
        (Exponent::Finite(_), Exponent::Finite(_)) => remquo_finite(rem, quo, lhs, rhs, rnd, stack),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use equator::assert;

    #[test]
    fn test_remainder_finite_direct_0() {
        let a = SmallFloat::from_parts(
            2,
            Sign::Pos,
            Exponent::Finite(10),
            utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let b = SmallFloat::from_parts(
            2,
            Sign::Pos,
            Exponent::Finite(1),
            utils::rev([0b1100_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );

        let mut q = [0];
        let mut c = SmallFloat::<1>::zero(3);
        remquo_finite_direct(
            &mut c,
            &mut q,
            &a,
            &b,
            Round::ToNearest,
            PodStack::new(bytemuck::cast_slice_mut(&mut [0u64; 100])),
        );
        assert!(c.to_f64(Round::ToNearest).0 == a.to_f64(Round::ToNearest).0 % b.to_f64(Round::ToNearest).0);
        assert!(q[0] as i64 == 341);
    }

    #[test]
    fn test_remainder_finite_direct_1() {
        let a = SmallFloat::from_parts(
            2,
            Sign::Pos,
            Exponent::Finite(10),
            utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let b = SmallFloat::from_parts(
            2,
            Sign::Neg,
            Exponent::Finite(1),
            utils::rev([0b1100_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );

        let mut q = [0];
        let mut c = SmallFloat::<1>::zero(3);
        remquo_finite_direct(
            &mut c,
            &mut q,
            &a,
            &b,
            Round::ToNearest,
            PodStack::new(bytemuck::cast_slice_mut(&mut [0u64; 100])),
        );
        assert!(c.to_f64(Round::ToNearest).0 == a.to_f64(Round::ToNearest).0 % b.to_f64(Round::ToNearest).0.abs());
        assert!(q[0] as i64 == 341);
    }

    #[test]
    fn test_remainder_finite_direct_2() {
        let a = SmallFloat::from_parts(
            2,
            Sign::Neg,
            Exponent::Finite(10),
            utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let b = SmallFloat::from_parts(
            2,
            Sign::Pos,
            Exponent::Finite(1),
            utils::rev([0b1100_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );

        let mut q = [0];
        let mut c = SmallFloat::<1>::zero(3);
        remquo_finite_direct(
            &mut c,
            &mut q,
            &a,
            &b,
            Round::ToNearest,
            PodStack::new(bytemuck::cast_slice_mut(&mut [0u64; 100])),
        );
        assert!(c.to_f64(Round::ToNearest).0 == a.to_f64(Round::ToNearest).0 % b.to_f64(Round::ToNearest).0 + b.to_f64(Round::ToNearest).0);
        assert!(q[0] as i64 == -342);
        assert!(a.to_f64(Round::ToNearest).0 - -342.0 * b.to_f64(Round::ToNearest).0 - c.to_f64(Round::ToNearest).0 == 0.0);
    }

    #[test]
    fn test_remainder_finite_direct_3() {
        let a = SmallFloat::from_parts(
            2,
            Sign::Neg,
            Exponent::Finite(10),
            utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let b = SmallFloat::from_parts(
            2,
            Sign::Neg,
            Exponent::Finite(1),
            utils::rev([0b1100_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );

        let mut q = [0];
        let mut c = SmallFloat::<1>::zero(3);
        remquo_finite_direct(
            &mut c,
            &mut q,
            &a,
            &b,
            Round::ToNearest,
            PodStack::new(bytemuck::cast_slice_mut(&mut [0u64; 100])),
        );
        assert!(
            c.to_f64(Round::ToNearest).0 == a.to_f64(Round::ToNearest).0 % b.to_f64(Round::ToNearest).0.abs() + b.to_f64(Round::ToNearest).0.abs()
        );
        assert!(q[0] as i64 == -342);
        assert!(a.to_f64(Round::ToNearest).0 - -342.0 * b.to_f64(Round::ToNearest).0.abs() - c.to_f64(Round::ToNearest).0 == 0.0);
    }

    #[test]
    fn test_remainder_finite_0() {
        let a = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(259),
            utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let b = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(1),
            utils::rev([0b1111_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );

        let mut q = [0];
        let mut c = SmallFloat::<1>::zero(3);

        remquo_finite(
            &mut c,
            &mut q,
            &a,
            &b,
            Round::ToNearest,
            PodStack::new(bytemuck::cast_slice_mut(&mut [0u64; 100])),
        );
        assert!(
            c.repr()
                == SmallFloat::from_parts(
                    3,
                    Sign::Pos,
                    Exponent::Finite(-1),
                    utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000])
                )
                .repr()
        );
        assert!(q == [2459565876494606882]);
    }

    #[test]
    fn test_remainder_finite_1() {
        let a = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(259),
            utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let b = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(1),
            utils::rev([0b1111_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );

        let mut q = [0, 0];
        let mut c = SmallFloat::<1>::zero(3);

        remquo_finite(
            &mut c,
            &mut q,
            &a,
            &b,
            Round::ToNearest,
            PodStack::new(bytemuck::cast_slice_mut(&mut [0u64; 100])),
        );
        assert!(
            c.repr()
                == SmallFloat::from_parts(
                    3,
                    Sign::Pos,
                    Exponent::Finite(-1),
                    utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000])
                )
                .repr()
        );
        assert!(
            q == utils::rev([
                0b10001000100010001000100010001000100010001000100010001000100010,
                0b0010001000100010001000100010001000100010001000100010001000100010
            ])
        );
    }

    #[test]
    fn test_remainder_finite_2() {
        let a = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(300),
            utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let b = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(1),
            utils::rev([0b1111_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );

        let mut q = [0, 0];
        let mut c = SmallFloat::<1>::zero(3);

        remquo_finite(
            &mut c,
            &mut q,
            &a,
            &b,
            Round::ToNearest,
            PodStack::new(bytemuck::cast_slice_mut(&mut [0u64; 100])),
        );
        assert!(
            c.repr()
                == SmallFloat::from_parts(
                    3,
                    Sign::Pos,
                    Exponent::Finite(0),
                    utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000])
                )
                .repr()
        );
        assert!(
            q == utils::rev([
                0b100010001000100010001000100010001000100010001000100010001000100,
                0b0100010001000100010001000100010001000100010001000100010001000100,
            ])
        );
    }

    #[test]
    fn test_remainder_finite_3() {
        let a = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(300),
            utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let b = SmallFloat::from_parts(
            4,
            Sign::Neg,
            Exponent::Finite(1),
            utils::rev([0b1111_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );

        let mut q = [0, 0];
        let mut c = SmallFloat::<1>::zero(3);

        remquo_finite(
            &mut c,
            &mut q,
            &a,
            &b,
            Round::ToNearest,
            PodStack::new(bytemuck::cast_slice_mut(&mut [0u64; 100])),
        );
        assert!(
            c.repr()
                == SmallFloat::from_parts(
                    3,
                    Sign::Pos,
                    Exponent::Finite(0),
                    utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000])
                )
                .repr()
        );
        assert!(
            q == utils::rev([
                0b100010001000100010001000100010001000100010001000100010001000100,
                0b0100010001000100010001000100010001000100010001000100010001000100,
            ])
        );
    }

    #[test]
    fn test_remainder_finite_4() {
        let a = SmallFloat::from_parts(
            4,
            Sign::Neg,
            Exponent::Finite(300),
            utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );
        let b = SmallFloat::from_parts(
            4,
            Sign::Pos,
            Exponent::Finite(1),
            utils::rev([0b1111_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000]),
        );

        let mut q = [0, 0];
        let mut c = SmallFloat::<1>::zero(4);

        remquo_finite(
            &mut c,
            &mut q,
            &a,
            &b,
            Round::ToNearest,
            PodStack::new(bytemuck::cast_slice_mut(&mut [0u64; 100])),
        );
        assert!(
            q == utils::rev([
                0b1011101110111011101110111011101110111011101110111011101110111011,
                0b1011101110111011101110111011101110111011101110111011101110111011,
            ])
        );
        assert!(
            c.repr()
                == SmallFloat::from_parts(
                    4,
                    Sign::Pos,
                    Exponent::Finite(1),
                    utils::rev([0b1011_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000])
                )
                .repr()
        );
    }
}
