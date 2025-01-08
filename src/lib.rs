use bytemuck::{Pod, Zeroable};
use core::{
    num::NonZeroU64,
    ops::{Shl, Shr},
};
use dyn_stack::{PodBuffer, PodStack, StackReq};
use equator::assert;
use reborrow::*;
#[allow(unused_imports)]
use utils::BinRepr as bin;

extern crate alloc;

mod copy;

mod add;
mod sub;

mod mul;

mod div;

mod sqrt;

mod remainder;

mod radix;

mod convert;

mod podstack;
pub use podstack::{temp_big_float_scratch, temp_big_float_uninit, temp_big_float_zero};

pub type Limb = u64;

pub mod utils {
    use crate::{BigFloat, Exponent, Limb, Sign};

    #[inline]
    pub const fn rev<const N: usize>(src: [Limb; N]) -> [Limb; N] {
        let mut dst = [0; N];
        let mut i = 0;
        while i < N {
            dst[i] = src[N - i - 1];
            i += 1;
        }
        dst
    }

    #[repr(transparent)]
    pub struct BinRepr(pub Limb);
    impl core::fmt::Debug for BinRepr {
        fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
            let width = crate::consts::LIMB_BITS as usize + 2;
            write!(f, "{:0>#0width$b}", self.0)
        }
    }

    #[repr(transparent)]
    pub struct FloatRepr {
        pub inner: BigFloat,
    }

    impl PartialEq for FloatRepr {
        fn eq(&self, other: &Self) -> bool {
            let prec_sign_exp =
                self.inner.precision_bits() == other.inner.precision_bits() && self.inner.sign_biased_exponent == other.inner.sign_biased_exponent;
            if !prec_sign_exp {
                false
            } else {
                match self.inner.exponent() {
                    Exponent::Zero | Exponent::NaN | Exponent::Inf => true,
                    Exponent::Finite(_) => self.inner.mantissa() == other.inner.mantissa(),
                }
            }
        }
    }

    impl core::fmt::Debug for FloatRepr {
        fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
            struct FloatRepr<'a> {
                precision_bits: u64,
                sign: Sign,
                exponent: Exponent,
                mantissa: &'a [BinRepr],
            }

            impl core::fmt::Debug for FloatRepr<'_> {
                fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
                    f.debug_struct("FloatRepr")
                        .field("precision_bits", &self.precision_bits)
                        .field("sign", &self.sign)
                        .field("exponent", &self.exponent)
                        .field("mantissa", &self.mantissa)
                        .finish()
                }
            }

            let mantissa = self.inner.mantissa();

            FloatRepr {
                precision_bits: self.inner.precision_bits(),
                sign: self.inner.sign(),
                exponent: self.inner.exponent(),
                mantissa: unsafe { core::slice::from_raw_parts(mantissa.as_ptr() as *const Limb as *const BinRepr, mantissa.len()) },
            }
            .fmt(f)
        }
    }
}

pub mod math {
    #[derive(Copy, Clone, Debug, PartialEq, Eq)]
    #[repr(i8)]
    pub enum Approx {
        Underflow = -2,
        LessThanExact = -1,
        Exact = 0,
        GreaterThanExact = 1,
        Overflow = 2,
    }

    #[derive(Copy, Clone, Debug, PartialEq, Eq)]
    pub enum Round {
        ToNearest,
        ToZero,
        AwayFromZero,
        Up,
        Down,
    }

    use super::*;

    /// TODO: docs
    pub use copy::abs;

    /// TODO: docs
    pub use copy::copy;

    /// TODO: docs
    #[inline]
    pub fn add(dst: &mut BigFloat, lhs: &BigFloat, rhs: &BigFloat, rnd: Round) -> Approx {
        if lhs.sign() == rhs.sign() {
            add::add_same_sign(dst, lhs.sign(), lhs, rhs, rnd)
        } else {
            sub::add_different_sign(dst, lhs.sign(), lhs, rhs, rnd)
        }
    }

    #[inline]
    pub fn sub(dst: &mut BigFloat, lhs: &BigFloat, rhs: &BigFloat, rnd: Round) -> Approx {
        if lhs.sign() == rhs.sign() {
            sub::add_different_sign(dst, lhs.sign(), lhs, rhs, rnd)
        } else {
            add::add_same_sign(dst, lhs.sign(), lhs, rhs, rnd)
        }
    }

    /// TODO: docs
    pub use mul::mul;

    /// TODO: docs
    pub use mul::mul_scratch;

    /// TODO: docs
    pub use mul::mul_add;

    /// TODO: docs
    pub use mul::mul_sub;

    /// TODO: docs
    pub use mul::negate_mul_add;

    /// TODO: docs
    pub use mul::negate_mul_sub;

    /// TODO: docs
    pub use div::div;

    /// TODO: docs
    pub use div::div_scratch;

    /// TODO: docs
    pub use sqrt::sqrt;

    /// TODO: docs
    pub use sqrt::sqrt_scratch;

    /// TODO: docs
    pub use remainder::remquo;

    /// TODO: docs
    pub use remainder::remquo_scratch;

    /// TODO: docs
    pub use convert::to_f64;

    /// TODO: docs
    pub use convert::from_f64;
}

#[derive(Copy, Clone, Debug)]
pub struct PrecisionCtx {
    precision_bits: NonZeroU64,
}

impl PrecisionCtx {
    #[track_caller]
    #[inline]
    pub fn new(precision_bits: u64) -> Self {
        assert!(all(
            precision_bits > 1,                         //
            precision_bits < BigFloat::max_precision(), //
        ));

        Self {
            precision_bits: NonZeroU64::new(precision_bits).unwrap(),
        }
    }

    #[inline]
    pub fn precision_bits(self) -> u64 {
        self.precision_bits.get()
    }

    pub fn add(self, lhs: &BigFloat, rhs: &BigFloat) -> BoxFloat {
        let mut out = BigFloat::zero(self.precision_bits());
        math::add(&mut out, lhs, rhs, Round::ToNearest);
        out
    }

    pub fn sub(self, lhs: &BigFloat, rhs: &BigFloat) -> BoxFloat {
        let mut out = BigFloat::zero(self.precision_bits());
        math::sub(&mut out, lhs, rhs, Round::ToNearest);
        out
    }

    pub fn abs(self, x: &BigFloat) -> BoxFloat {
        let mut out = BigFloat::zero(self.precision_bits());
        math::abs(&mut out, x, Round::ToNearest);
        out
    }

    pub fn mul(self, lhs: &BigFloat, rhs: &BigFloat) -> BoxFloat {
        let mut out = BigFloat::zero(self.precision_bits());
        math::mul(
            &mut out,
            lhs,
            rhs,
            Round::ToNearest,
            PodStack::new(&mut PodBuffer::new(math::mul_scratch(
                self.precision_bits(),
                lhs.precision_bits(),
                rhs.precision_bits(),
            ))),
        );
        out
    }

    pub fn div(self, lhs: &BigFloat, rhs: &BigFloat) -> BoxFloat {
        let mut out = BigFloat::zero(self.precision_bits());
        math::div(
            &mut out,
            lhs,
            rhs,
            Round::ToNearest,
            PodStack::new(&mut PodBuffer::new(math::div_scratch(
                self.precision_bits(),
                lhs.precision_bits(),
                rhs.precision_bits(),
            ))),
        );
        out
    }

    pub fn sqrt(self, x: &BigFloat) -> BoxFloat {
        let mut out = BigFloat::zero(self.precision_bits());
        math::sqrt(
            &mut out,
            x,
            Round::ToNearest,
            PodStack::new(&mut PodBuffer::new(math::sqrt_scratch(self.precision_bits(), x.precision_bits()))),
        );
        out
    }

    pub fn rem(self, lhs: &BigFloat, rhs: &BigFloat) -> BoxFloat {
        self.remquo(&mut [], lhs, rhs)
    }

    pub fn remquo(self, quo: &mut [Limb], lhs: &BigFloat, rhs: &BigFloat) -> BoxFloat {
        let mut out = BigFloat::zero(self.precision_bits());
        let quo_len = quo.len();
        math::remquo(
            &mut out,
            quo,
            lhs,
            rhs,
            Round::ToNearest,
            PodStack::new(&mut PodBuffer::new(math::remquo_scratch(
                self.precision_bits(),
                quo_len,
                lhs.precision_bits(),
                rhs.precision_bits(),
            ))),
        );
        out
    }
}

#[repr(C)]
pub struct BigFloat {
    pub sign_biased_exponent: u64,
    __precision_bits: Option<NonZeroU64>,
    __mantissa: [Limb],
}

#[repr(C)]
#[derive(Copy, Clone)]
pub struct SmallFloat<const N: usize> {
    pub sign_biased_exponent: u64,
    __precision_bits: Option<NonZeroU64>,
    __mantissa: [Limb; N],
}

unsafe impl<const N: usize> Zeroable for SmallFloat<N> {}
unsafe impl<const N: usize> Pod for SmallFloat<N> {}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Sign {
    Neg,
    Pos,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Exponent {
    Zero,
    Finite(i64),
    Inf,
    NaN,
}

impl Round {
    pub const ALL: &'static [Self] = &[Self::ToNearest, Self::ToZero, Self::AwayFromZero, Self::Up, Self::Down];
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

use math::{Approx, Round};

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
    #[must_use]
    pub const fn neg(self) -> Self {
        match self {
            Sign::Neg => Self::Pos,
            Sign::Pos => Self::Neg,
        }
    }
}

pub struct BoxFloat {
    inner: alloc::boxed::Box<BigFloat>,
}

impl core::ops::Deref for BoxFloat {
    type Target = BigFloat;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &*self.inner
    }
}

impl core::ops::DerefMut for BoxFloat {
    #[inline]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut *self.inner
    }
}

impl Clone for BoxFloat {
    fn clone(&self) -> Self {
        let mut tmp = BigFloat::zero(self.inner.precision_bits());
        tmp.inner.sign_biased_exponent = self.inner.sign_biased_exponent;
        tmp.inner.__mantissa.copy_from_slice(&self.inner.__mantissa);
        tmp
    }
}

pub mod consts {
    use crate::Limb;

    pub const MAX_EXPONENT_INCLUSIVE: i64 = i64::MAX / 2 - 1;
    pub const MIN_EXPONENT_INCLUSIVE: i64 = -(i64::MAX / 2 - 1);

    pub const EXPONENT_BIAS: i64 = i64::MAX / 2;
    pub const BIASED_EXPONENT_NAN: u64 = i64::MAX as u64;
    pub const BIASED_EXPONENT_INF: u64 = i64::MAX as u64 - 1;

    pub const SIGN_SHIFT: u64 = 63;
    pub const SIGN_BIT: u64 = 1u64 << SIGN_SHIFT;
    pub const BIASED_EXPONENT_MASK: u64 = (1u64 << 63) - 1;

    pub const LIMB_ZERO: Limb = 0;
    pub const LIMB_ONE: Limb = 1;
    pub const LIMB_BITS: u64 = Limb::BITS as u64;
    pub const LIMB_HIGH_BIT: Limb = 0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000;
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

impl BigFloat {
    #[must_use]
    #[track_caller]
    #[inline]
    pub fn zero(precision_bits: u64) -> BoxFloat {
        use alloc::alloc::*;

        assert!(all(precision_bits > 1, precision_bits < u64::MAX - Limb::BITS as u64,));

        let precision_bits_usize: usize = precision_bits.try_into().unwrap();

        let nlimbs = precision_bits_usize.div_ceil(Limb::BITS as usize);
        let layout = Layout::array::<Limb>(nlimbs)
            .and_then(|tail| Layout::new::<u64>().extend(Layout::new::<u64>()).unwrap().0.extend(tail))
            .unwrap()
            .0
            .pad_to_align();
        let ptr = unsafe { alloc_zeroed(layout) };
        if ptr.is_null() {
            handle_alloc_error(layout);
        }

        let ptr = core::ptr::slice_from_raw_parts_mut(ptr, nlimbs) as *mut BigFloat;
        unsafe { (*ptr).__precision_bits = Some(NonZeroU64::new_unchecked(precision_bits)) };
        BoxFloat {
            inner: unsafe { alloc::boxed::Box::from_raw(ptr) },
        }
    }

    #[inline]
    pub const fn precision_bits(&self) -> u64 {
        match self.__precision_bits {
            Some(precision_bits) => precision_bits.get(),
            None => self.__mantissa.len() as u64 * Limb::BITS as u64,
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
    pub const fn mantissa_len(&self) -> usize {
        ((self.precision_bits() + (Limb::BITS as u64 - 1)) / Limb::BITS as u64) as usize
    }

    #[inline]
    pub const fn mantissa(&self) -> &[Limb] {
        self.__mantissa.split_at(self.mantissa_len()).0
    }

    #[inline]
    pub fn mantissa_mut(&mut self) -> &mut [Limb] {
        self.__mantissa.split_at_mut(self.mantissa_len()).0
    }

    #[inline]
    pub const fn full_mantissa(&self) -> &[Limb] {
        &self.__mantissa
    }

    #[inline]
    pub fn full_mantissa_mut(&mut self) -> &mut [Limb] {
        &mut self.__mantissa
    }

    #[inline]
    pub fn copy_from(&mut self, src: &BigFloat, rnd: Round) -> Approx {
        copy::copy(self, src, rnd)
    }

    #[inline]
    pub fn to_f64(&self, rnd: Round) -> (f64, Approx) {
        convert::to_f64(self, rnd)
    }

    #[inline]
    pub fn repr(&self) -> &utils::FloatRepr {
        unsafe { &*(self as *const BigFloat as *const utils::FloatRepr) }
    }

    #[inline]
    pub const fn max_precision() -> u64 {
        u64::MAX - consts::LIMB_BITS
    }
}

impl<const N: usize> SmallFloat<N> {
    #[inline]
    #[track_caller]
    pub const fn max_precision() -> u64 {
        if N == 0 {
            panic!();
        }

        let max = u64::MAX - consts::LIMB_BITS;
        match (N as u64).checked_mul(consts::LIMB_BITS) {
            Some(p) => {
                if p < max {
                    p
                } else {
                    max
                }
            }
            None => panic!(),
        }
    }

    #[inline]
    #[track_caller]
    pub const fn from_parts(precision_bits: u64, sign: Sign, exponent: Exponent, mantissa: [Limb; N]) -> Self {
        if N == 0 {
            panic!();
        }
        if precision_bits <= 1 || precision_bits >= Self::max_precision() {
            panic!()
        }

        Self {
            sign_biased_exponent: make_sign_and_biased_exponent(sign, exponent),
            __precision_bits: Some(unsafe { NonZeroU64::new_unchecked(precision_bits) }),
            __mantissa: mantissa,
        }
    }

    #[inline]
    #[track_caller]
    pub const fn zero(precision_bits: u64) -> Self {
        if precision_bits == 0 || precision_bits >= Self::max_precision() {
            panic!()
        }

        Self {
            sign_biased_exponent: 0,
            __precision_bits: Some(unsafe { NonZeroU64::new_unchecked(precision_bits) }),
            __mantissa: [0; N],
        }
    }

    #[inline]
    #[must_use]
    pub const fn as_ref(&self) -> &BigFloat {
        unsafe { &*(core::ptr::slice_from_raw_parts(self as *const SmallFloat<N> as *const Limb, N) as *const BigFloat) }
    }

    #[inline]
    #[must_use]
    pub fn as_mut(&mut self) -> &mut BigFloat {
        unsafe { &mut *(core::ptr::slice_from_raw_parts_mut(self as *mut SmallFloat<N> as *mut Limb, N) as *mut BigFloat) }
    }
}

impl<const N: usize> core::ops::Deref for SmallFloat<N> {
    type Target = BigFloat;

    #[inline]
    #[track_caller]
    fn deref(&self) -> &Self::Target {
        self.as_ref()
    }
}

impl<const N: usize> core::ops::DerefMut for SmallFloat<N> {
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
        let x = BigFloat::zero(2);
        _ = x.clone();
        let a = SmallFloat::from_parts(3, Sign::Pos, Exponent::Finite(3), [5, 9u64]);
        let mut b = SmallFloat::from_parts(3, Sign::Pos, Exponent::Finite(3), [192, 28u64]);

        _ = a.as_ref();
        _ = b.as_mut();
    }
}
