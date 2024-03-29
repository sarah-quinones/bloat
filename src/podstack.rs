use super::*;
use dyn_stack::{PodStack, SizeOverflow, StackReq};
use equator::assert;

#[inline]
pub fn temp_big_float_req(precision_bits: u64) -> Result<StackReq, SizeOverflow> {
    let nlimbs = (precision_bits as usize).div_ceil(Limb::BITS as usize);
    StackReq::try_new::<u64>(2 + nlimbs)
}

#[inline]
#[track_caller]
pub fn temp_big_float_uninit(precision_bits: u64, stack: PodStack<'_>) -> (&mut BigFloat, PodStack<'_>) {
    assert!(core::mem::size_of::<Limb>() == 8);
    assert!(all(precision_bits > 1, precision_bits < u64::MAX - Limb::BITS as u64,));
    let nlimbs = (precision_bits as usize).div_ceil(Limb::BITS as usize);
    let (buf, stack) = stack.make_raw::<u64>(2 + nlimbs);
    let buf = buf.as_mut_ptr();
    let buf = core::ptr::slice_from_raw_parts_mut(buf, nlimbs) as *mut BigFloat;
    let buf = unsafe { &mut *buf };
    buf.__precision_bits = Some(unsafe { NonZeroU64::new_unchecked(precision_bits) });
    (buf, stack)
}

#[inline]
#[track_caller]
pub fn temp_big_float_zero(precision_bits: u64, stack: PodStack<'_>) -> (&mut BigFloat, PodStack<'_>) {
    let (buf, stack) = temp_big_float_uninit(precision_bits, stack);
    buf.sign_biased_exponent = 0;
    buf.mantissa_mut().fill(consts::LIMB_ZERO);
    (buf, stack)
}
