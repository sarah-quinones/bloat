use std::num::NonZeroU64;

use equator::assert;

extern crate alloc;

pub mod utils {
    #[inline]
    pub const fn rev<L: crate::Limb, const N: usize>(src: [L; N]) -> [L; N] {
        let mut dst = [L::ZERO; N];
        let mut i = 0;
        while i < N {
            dst[i] = src[N - i - 1];
            i += 1;
        }
        dst
    }
}

mod add;
mod copy;

pub mod math {
    pub use crate::copy::{abs, copy};
}

pub trait Limb:
    Copy
    + core::cmp::PartialEq
    + core::cmp::PartialOrd
    + core::ops::Not<Output = Self>
    + core::ops::BitAnd<Self, Output = Self>
    + core::ops::BitOr<Self, Output = Self>
    + core::ops::BitXor<Self, Output = Self>
    + core::ops::BitAndAssign<Self>
    + core::ops::BitOrAssign<Self>
    + core::ops::BitXorAssign<Self>
{
    const BITS: u64;
    const ZERO: Self;
    const ONE: Self;

    fn from_bit(bit: bool) -> Self;

    fn overflowing_add(self, rhs: Self) -> (Self, bool);
    fn overflowing_sub(self, rhs: Self) -> (Self, bool);

    fn wrapping_add(self, rhs: Self) -> Self;
    fn wrapping_sub(self, rhs: Self) -> Self;
    fn widening_mul(self, rhs: Self) -> (Self, Self);
    fn narrowing_div(lhs: (Self, Self), rhs: Self) -> Self;
    fn shl(self, amount: u64) -> Self;
    fn shr(self, amount: u64) -> Self;
}

impl Limb for u32 {
    const BITS: u64 = Self::BITS as u64;
    const ZERO: Self = 0;
    const ONE: Self = 1;

    #[inline]
    fn from_bit(bit: bool) -> Self {
        bit as Self
    }

    #[inline]
    fn overflowing_add(self, rhs: Self) -> (Self, bool) {
        Self::overflowing_add(self, rhs)
    }
    #[inline]
    fn overflowing_sub(self, rhs: Self) -> (Self, bool) {
        Self::overflowing_sub(self, rhs)
    }

    #[inline]
    fn wrapping_add(self, rhs: Self) -> Self {
        Self::wrapping_add(self, rhs)
    }

    #[inline]
    fn wrapping_sub(self, rhs: Self) -> Self {
        Self::wrapping_sub(self, rhs)
    }

    #[inline]
    fn widening_mul(self, rhs: Self) -> (Self, Self) {
        let wide = self as u64 * rhs as u64;
        (wide as u32, (wide >> 32) as u32)
    }

    #[inline]
    #[track_caller]
    fn narrowing_div(lhs: (Self, Self), rhs: Self) -> Self {
        let wide = lhs.0 as u64 + ((lhs.1 as u64) << 32);
        (wide / rhs as u64) as u32
    }

    #[inline]
    fn shl(self, amount: u64) -> Self {
        self << amount
    }

    #[inline]
    fn shr(self, amount: u64) -> Self {
        self >> amount
    }
}

impl Limb for u64 {
    const BITS: u64 = Self::BITS as u64;
    const ZERO: Self = 0;
    const ONE: Self = 1;

    #[inline]
    fn from_bit(bit: bool) -> Self {
        bit as Self
    }

    #[inline]
    fn overflowing_add(self, rhs: Self) -> (Self, bool) {
        Self::overflowing_add(self, rhs)
    }
    #[inline]
    fn overflowing_sub(self, rhs: Self) -> (Self, bool) {
        Self::overflowing_sub(self, rhs)
    }

    #[inline]
    fn wrapping_add(self, rhs: Self) -> Self {
        Self::wrapping_add(self, rhs)
    }

    #[inline]
    fn wrapping_sub(self, rhs: Self) -> Self {
        Self::wrapping_sub(self, rhs)
    }

    #[inline]
    fn widening_mul(self, rhs: Self) -> (Self, Self) {
        let wide = self as u128 * rhs as u128;
        (wide as u64, (wide >> 64) as u64)
    }

    #[inline]
    #[track_caller]
    fn narrowing_div(lhs: (Self, Self), rhs: Self) -> Self {
        let wide = lhs.0 as u128 + ((lhs.1 as u128) << 64);
        (wide / rhs as u128) as u64
    }

    #[inline]
    fn shl(self, amount: u64) -> Self {
        self << amount
    }

    #[inline]
    fn shr(self, amount: u64) -> Self {
        self >> amount
    }
}

#[repr(C)]
pub struct BigFloat<L> {
    sign_biased_exponent: u64,
    __precision_bits: Option<NonZeroU64>,
    mantissa: [L],
}

#[repr(C)]
#[derive(Debug)]
pub struct SmallFloat<L, const N: usize> {
    sign_biased_exponent: u64,
    __precision_bits: Option<NonZeroU64>,
    mantissa: [L; N],
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Sign {
    Neg,
    Pos,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Exponent {
    Finite(i64),
    Zero,
    NaN,
    Inf,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Round {
    ToNearest,
    ToZero,
    AwayFromZero,
    Up,
    Down,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum RoundKnownSign {
    ToNearest,
    ToZero,
    AwayFromZero,
}

impl Round {
    #[inline]
    fn with_sign(self, sign: Sign) -> RoundKnownSign {
        match self {
            Round::ToNearest => RoundKnownSign::ToNearest,
            Round::ToZero => RoundKnownSign::ToZero,
            Round::AwayFromZero => RoundKnownSign::AwayFromZero,
            Round::Up => {
                if sign.is_negative() {
                    RoundKnownSign::ToZero
                } else {
                    RoundKnownSign::AwayFromZero
                }
            }
            Round::Down => {
                if sign.is_positive() {
                    RoundKnownSign::ToZero
                } else {
                    RoundKnownSign::AwayFromZero
                }
            }
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(i8)]
pub enum Approx {
    LessThanExact = -1,
    Exact = 0,
    GreaterThanExact = 1,
}

impl Approx {
    #[inline]
    fn from_sign(sign: Sign) -> Self {
        if sign.is_negative() {
            Self::LessThanExact
        } else {
            Self::GreaterThanExact
        }
    }
}

impl Sign {
    #[inline]
    pub const fn is_positive(self) -> bool {
        matches!(self, Self::Pos)
    }
    #[inline]
    pub const fn is_negative(self) -> bool {
        matches!(self, Self::Neg)
    }

    #[inline]
    pub const fn neg(self) -> Self {
        match self {
            Sign::Neg => Self::Pos,
            Sign::Pos => Self::Neg,
        }
    }
}

pub type BoxFloat<T> = alloc::boxed::Box<BigFloat<T>>;

impl<L: Limb> Clone for BoxFloat<L> {
    fn clone(&self) -> Self {
        let mut tmp = BigFloat::zero(self.precision_bits());
        tmp.sign_biased_exponent = self.sign_biased_exponent;
        tmp.mantissa.copy_from_slice(&self.mantissa);
        tmp
    }
}

pub mod consts {
    pub const MAX_EXPONENT_INCLUSIVE: i64 = i64::MAX / 2 - 1;
    pub const MIN_EXPONENT_INCLUSIVE: i64 = -(i64::MAX / 2 - 1);

    pub const EXPONENT_BIAS: i64 = i64::MAX / 2;
    pub const BIASED_EXPONENT_NAN: u64 = i64::MAX as u64;
    pub const BIASED_EXPONENT_INF: u64 = i64::MAX as u64 - 1;

    pub const SIGN_SHIFT: u64 = 63;
    pub const SIGN_BIT: u64 = 1u64 << SIGN_SHIFT;
    pub const BIASED_EXPONENT_MASK: u64 = (1u64 << 63) - 1;
}

#[inline]
pub const fn make_sign_and_biased_exponent(sign: Sign, exponent: Exponent) -> u64 {
    let sign = (sign.is_negative() as u64) << 63;
    let exp = match exponent {
        Exponent::Finite(exp) => {
            if exp > consts::MAX_EXPONENT_INCLUSIVE || exp < consts::MIN_EXPONENT_INCLUSIVE {
                panic!()
            }
            exp.wrapping_add(consts::EXPONENT_BIAS) as u64
        }
        Exponent::Zero => 0,
        Exponent::NaN => consts::BIASED_EXPONENT_NAN,
        Exponent::Inf => consts::BIASED_EXPONENT_INF,
    };
    sign | exp
}

impl<L: Limb> BigFloat<L> {
    #[must_use]
    #[track_caller]
    #[inline]
    pub fn zero(precision_bits: u64) -> BoxFloat<L> {
        use alloc::alloc::*;

        assert!(precision_bits > 0);

        let precision_bits_usize: usize = precision_bits.try_into().unwrap();

        let nlimbs = precision_bits_usize.div_ceil(L::BITS as usize);
        let layout = Layout::array::<L>(nlimbs)
            .and_then(|tail| {
                Layout::new::<u64>()
                    .extend(Layout::new::<u64>())
                    .unwrap()
                    .0
                    .extend(tail)
            })
            .unwrap()
            .0
            .pad_to_align();
        let ptr = unsafe { alloc_zeroed(layout) };
        if ptr.is_null() {
            handle_alloc_error(layout);
        }

        let ptr = core::ptr::slice_from_raw_parts_mut(ptr, nlimbs) as *mut BigFloat<L>;
        unsafe { (*ptr).__precision_bits = Some(NonZeroU64::new_unchecked(precision_bits)) };
        unsafe { alloc::boxed::Box::from_raw(ptr) }
    }

    #[inline]
    pub const fn precision_bits(&self) -> u64 {
        match self.__precision_bits {
            Some(precision_bits) => precision_bits.get(),
            None => self.mantissa.len() as u64 * L::BITS,
        }
    }

    #[inline]
    pub const fn sign(&self) -> Sign {
        if (self.sign_biased_exponent >> 63) == 1 {
            Sign::Neg
        } else {
            Sign::Pos
        }
    }

    #[inline]
    pub const fn exponent(&self) -> Exponent {
        match self.sign_biased_exponent & ((1 << 63) - 1) {
            0 => Exponent::Zero,
            consts::BIASED_EXPONENT_NAN => Exponent::NaN,
            consts::BIASED_EXPONENT_INF => Exponent::Inf,
            biased => Exponent::Finite((biased as i64).wrapping_sub(consts::EXPONENT_BIAS)),
        }
    }

    #[inline]
    pub const fn mantissa(&self) -> &[L] {
        &self.mantissa
    }

    #[inline]
    pub fn copy_from(&mut self, src: &BigFloat<L>, rnd: Round) -> Approx {
        copy::copy(self, src, rnd)
    }
}

impl<L: Limb, const N: usize> SmallFloat<L, N> {
    #[inline]
    #[track_caller]
    pub const fn from_parts(
        sign: Sign,
        exponent: Exponent,
        precision_bits: u64,
        mantissa: [L; N],
    ) -> Self {
        if precision_bits == 0 || precision_bits > N as u64 * L::BITS {
            panic!()
        }

        Self {
            sign_biased_exponent: make_sign_and_biased_exponent(sign, exponent),
            __precision_bits: Some(unsafe { NonZeroU64::new_unchecked(precision_bits) }),
            mantissa,
        }
    }

    #[inline]
    #[track_caller]
    pub const fn zero(precision_bits: u64) -> Self {
        if precision_bits == 0 || precision_bits > N as u64 * L::BITS {
            panic!()
        }

        Self {
            sign_biased_exponent: 0,
            __precision_bits: Some(unsafe { NonZeroU64::new_unchecked(precision_bits) }),
            mantissa: [L::ZERO; N],
        }
    }

    #[inline]
    #[track_caller]
    #[must_use]
    pub const fn as_ref(&self) -> &BigFloat<L> {
        if N == 0 {
            panic!();
        }
        unsafe {
            &*(core::ptr::slice_from_raw_parts(self as *const SmallFloat<L, N> as *const L, N)
                as *const BigFloat<L>)
        }
    }

    #[inline]
    #[track_caller]
    #[must_use]
    pub fn as_mut(&mut self) -> &mut BigFloat<L> {
        if N == 0 {
            panic!();
        }
        unsafe {
            &mut *(core::ptr::slice_from_raw_parts_mut(self as *mut SmallFloat<L, N> as *mut L, N)
                as *mut BigFloat<L>)
        }
    }
}

impl<L: Limb, const N: usize> core::ops::Deref for SmallFloat<L, N> {
    type Target = BigFloat<L>;

    #[inline]
    #[track_caller]
    fn deref(&self) -> &Self::Target {
        self.as_ref()
    }
}

impl<L: Limb, const N: usize> core::ops::DerefMut for SmallFloat<L, N> {
    #[inline]
    #[track_caller]
    fn deref_mut(&mut self) -> &mut Self::Target {
        self.as_mut()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_unsafe() {
        let x = BigFloat::<u64>::zero(2);
        _ = x.clone();
        let a = SmallFloat::from_parts(Sign::Pos, Exponent::Finite(3), 3, [5, 9u64]);
        let mut b = SmallFloat::from_parts(Sign::Pos, Exponent::Finite(3), 3, [192, 28u64]);

        _ = a.as_ref();
        _ = b.as_mut();
    }
}
