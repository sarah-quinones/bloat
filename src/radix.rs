use super::*;
use core::fmt;

pub fn to_decimal_scratch(x_bits: u64) -> StackReq {
    _ = x_bits;
    StackReq::EMPTY
}

pub fn to_decimal(f: &mut dyn fmt::Write, x: &BigFloat, stack: &mut PodStack) -> fmt::Result {
    _ = stack;
    write!(f, "{}", convert::to_f64(x, Round::ToNearest).0)
}

impl fmt::Display for BigFloat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        to_decimal(
            f,
            self,
            PodStack::new(&mut dyn_stack::PodBuffer::new(to_decimal_scratch(self.precision_bits()))),
        )
    }
}

impl<const N: usize> fmt::Display for SmallFloat<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        (**self).fmt(f)
    }
}
