use super::*;
use core::fmt;

pub fn to_decimal_req(x_bits: u64) -> Result<StackReq, SizeOverflow> {
    _ = x_bits;
    Ok(StackReq::empty())
}

pub fn to_decimal(f: &mut dyn fmt::Write, x: &BigFloat, stack: PodStack<'_>) -> fmt::Result {
    _ = stack;
    write!(f, "{}", convert::to_f64(x, Round::ToNearest).0)
}

impl fmt::Display for BigFloat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        to_decimal(
            f,
            self,
            PodStack::new(&mut dyn_stack::GlobalPodBuffer::new(to_decimal_req(self.precision_bits()).unwrap())),
        )
    }
}

impl<const N: usize> fmt::Display for SmallFloat<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        (**self).fmt(f)
    }
}
