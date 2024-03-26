// use super::*;
// use equator::debug_assert;

// pub fn add_same_sign(dst: &mut BigFloat, lhs: &BigFloat, rhs: &BigFloat, rnd: Round) -> Approx {
//     debug_assert!(lhs.sign_is_negative() == rhs.sign_is_negative());

//     let (lhs, rhs) = if (lhs.sign_biased_exponent & (i64::MAX as u64))
//         > (rhs.sign_biased_exponent & (i64::MAX as u64))
//     {
//         (lhs, rhs)
//     } else {
//         (rhs, lhs)
//     };

//     let lhs_exp = lhs.exponent();
//     let rhs_exp = rhs.exponent();
//     if rhs_exp == Exponent::Zero || lhs_exp == Exponent::NaN || lhs_exp == Exponent::Inf {
//         return dst.copy_from(lhs, rnd);
//     }

//     let is_negative = lhs.sign_is_negative();
//     let (Exponent::FiniteNonZero(lhs_exp), Exponent::FiniteNonZero(rhs_exp)) = (lhs_exp, rhs_exp)
//     else {
//         unreachable!()
//     };

//     let shift = (lhs_exp as u64).wrapping_sub(rhs_exp as u64);
//     let lhs_len = lhs.mantissa.len();
//     let rhs_len = rhs.mantissa.len();

//     todo!()
// }
