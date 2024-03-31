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
