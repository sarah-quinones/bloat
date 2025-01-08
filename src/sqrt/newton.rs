use super::*;
use equator::assert;

pub fn isqrt_scratch(n_len: usize, lshift: u64) -> StackReq {
    let nbits = n_len as u64 * (consts::LIMB_BITS / 2) + lshift.div_ceil(2);
    let nlimbs = nbits / consts::LIMB_BITS;
    let out_len = nlimbs as usize;

    // let div_lhs_len = n_len + 1 + lshift.div_ceil(consts::LIMB_BITS) as usize;
    StackReq::all_of(&[
        StackReq::new::<Limb>(n_len),
        StackReq::new::<Limb>(1 + lshift.div_ceil(consts::LIMB_BITS) as usize),
        StackReq::new::<Limb>(out_len),
    ])
}

pub fn isqrt(out: &mut [Limb], n: &[Limb], lshift: u64, stack: &mut PodStack) {
    let (div_lhs, stack) = stack.make_raw::<Limb>(n.len() + 1 + lshift.div_ceil(consts::LIMB_BITS) as usize);
    let (div_quo, mut stack) = stack.make_raw::<Limb>(out.len());

    assert!(all(
        consts::LIMB_BITS % 2 == 0,            //
    ));
    let nbits = n.len() as u64 * (consts::LIMB_BITS / 2) + lshift.div_ceil(2);
    let nlimbs = nbits / consts::LIMB_BITS;

    assert!(all(
        nbits % consts::LIMB_BITS == 0, //
        out.len() == nlimbs as usize,   //
    ));

    out.fill(!consts::LIMB_ZERO);
    loop {
        // need to shift by e.g. 75
        // my "limbs" (groups of bits) contain 64 bits
        // shift by 64, then shift by 75 - 64 = 11

        // we call shift by 64 a large shift, because we shift by an entire limb
        // 11 is a small shift, because it's within a single limb

        //       [n1 n0]
        // [m3 m2 m1 m0]
        // m0 = 0
        // m1 = (n0 << 11)
        // m2 = (n1 << 11) | (n0 >> 53)
        // m3 =              (n1 >> 53)

        div_lhs.fill(consts::LIMB_ZERO);
        let large_lshift = (lshift / consts::LIMB_BITS) as usize;
        let small_lshift = lshift % consts::LIMB_BITS;

        div_lhs[large_lshift] = n[0].shl(small_lshift);
        if small_lshift == 0 {
            for i in 1..n.len() {
                div_lhs[large_lshift + i] = n[i].shl(small_lshift)
            }
        } else {
            for i in 1..n.len() {
                div_lhs[large_lshift + i] = n[i].shl(small_lshift) | n[i - 1].shr(consts::LIMB_BITS - small_lshift)
            }
            div_lhs[large_lshift + n.len()] = n[n.len() - 1].shr(consts::LIMB_BITS - small_lshift);
        }

        // compute (n // x) in div_quo
        div::idivrem_normalized(div_quo, None, div_lhs, &*out, stack.rb_mut());

        // compute (x + n // x) // 2 in x
        let mut carry;
        let mut old = out[0];
        let mut same_as_old = true;

        (out[0], carry) = out[0].overflowing_add(div_quo[0]);
        out[0] >>= 1;
        for i in 1..out.len() {
            let (sum, carry0) = out[i].overflowing_add(div_quo[i]);
            let (sum, carry1) = sum.overflowing_add(carry as Limb);

            carry = carry0 | carry1;
            out[i - 1] |= sum << (consts::LIMB_BITS - 1);
            same_as_old &= out[i - 1] == old;
            old = out[i];
            out[i] = sum >> 1;
        }
        let len = out.len();
        out[len - 1] |= (carry as Limb) << (consts::LIMB_BITS - 1);
        same_as_old &= out[len - 1] == old;

        if same_as_old {
            break;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use equator::assert;

    #[test]
    fn test_isqrt_0() {
        let n = utils::rev([0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000, 0]);
        let mut out = [0];

        isqrt(&mut out, &n, 0, PodStack::new(bytemuck::cast_slice_mut(&mut [0u64; 100])));
        assert!(out == [13043817825332782212]);
    }

    #[test]
    fn test_isqrt_1() {
        let n = utils::rev([
            0b1001010110010100111010101001110100111111111111101100000110011101,
            0b1100001101101110000001100000010110100011111100110110111111000011,
            0b0000110101011100011100011111011011000101000011110100111000101011,
            0b1100111010100101011001111101111110000111010001100111011111111101,
        ]);

        let mut out = [0, 0];
        isqrt(&mut out, &n, 0, PodStack::new(bytemuck::cast_slice_mut(&mut [0u64; 100])));
        assert!(
            out == utils::rev([
                0b1100001110101111100011011110111011011101000010001110011111000100,
                0b0111111101100011110101000100010001001000111010111101101001001011,
            ])
        );
    }

    #[test]
    fn test_isqrt_2() {
        let n = utils::rev([
            0b1001010110010100111010101001110100111111111111101100000110011101,
            0b1100001101101110000001100000010110100011111100110110111111000011,
            0b0000110101011100011100011111011011000101000011110100111000101011,
            0b1100111010100101011001111101111110000111010001100111011111111101,
        ]);
        let mut out = [0, 0, 0];
        isqrt(&mut out, &n, 127, PodStack::new(bytemuck::cast_slice_mut(&mut [0u64; 100])));
        assert!(
            out == utils::rev([
                0b1000101001011110111001111111111101111000101101011000001010000111,
                0b0011001001101001000101101010111111000100010001100100101110000100,
                0b0001001101111100010011001000010001101100010111100100010011101100,
            ])
        );
    }
}
