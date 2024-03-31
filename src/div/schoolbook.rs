use super::*;
use equator::debug_assert;

fn widening_mul_acc(lhs: Limb, rhs: Limb, acc: Limb, carry: Limb) -> (Limb, Limb) {
    let wide = acc as u128 + carry as u128 + lhs as u128 * rhs as u128;
    (wide as u64, (wide >> 64) as u64)
}

fn to_wide(lo: u64, hi: u64) -> u128 {
    lo as u128 | ((hi as u128) << 64)
}

#[inline]
pub fn divrem_bigint(quo: &mut [Limb], rem: Option<&mut [Limb]>, lhs: &mut [Limb], rhs: &[Limb]) {
    debug_assert!(rhs[rhs.len() - 1] >> (consts::LIMB_BITS - 1) == 1);
    debug_assert!(lhs[lhs.len() - 1] == 0);

    let m = quo.len() - 1;
    let n = rhs.len();

    let u = lhs;
    let v = rhs;

    quo.fill(consts::LIMB_ZERO);

    let mut j = m + 1;
    while j != 0 {
        j -= 1;

        let u_wide = to_wide(u[n + j - 1], u[n + j]);
        let mut qhat = u_wide / v[n - 1] as u128;
        let mut rhat = u_wide % v[n - 1] as u128;
        loop {
            if qhat == to_wide(0, 1) || (n > 1 && qhat * v[n - 2] as u128 > to_wide(u[n + j - 2], rhat as u64)) {
                qhat -= 1;
                rhat += v[n - 1] as u128;
                if rhat >= to_wide(0, 1) {
                    break;
                }
            } else {
                break;
            }
        }
        let qhat = qhat as u64;

        let qhat_c = qhat.wrapping_neg();
        let mut borrow = false;
        let mut carry = false;

        if qhat != 0 {
            let mut carry_u64 = 0u64;
            for i in 0..n {
                (u[i + j], carry_u64) = widening_mul_acc(qhat_c, v[i], u[i + j], carry_u64);
            }
            (u[n + j], carry) = u[n + j].overflowing_add(carry_u64);

            for i in 1..n + 1 {
                let (diff, borrow0) = u[i + j].overflowing_sub(v[i - 1]);
                let (diff, borrow1) = diff.overflowing_sub(borrow as u64);
                u[i + j] = diff;
                borrow = borrow0 | borrow1;
            }
        }
        match (carry, borrow) {
            (true, true) => borrow = false,
            (false, false) | (false, true) => {}
            (true, false) => panic!(),
        }

        if borrow {
            if qhat != 1 {
                quo[j] = qhat - 1;
            }
            let mut carry = false;
            for i in 0..n {
                let (sum, carry0) = u[i + j].overflowing_sub(v[i]);
                let (sum, carry1) = sum.overflowing_add(carry as u64);
                u[i + j] = sum;
                carry = carry0 | carry1;
            }
            u[n + j] = u[n + j].wrapping_add(carry as u64);
        } else {
            if qhat != 0 {
                quo[j] = qhat;
            }
        }
    }

    if let Some(rem) = rem {
        rem.copy_from_slice(&u[..rem.len()]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use equator::assert;

    #[test]
    fn test_div_bigint_0() {
        let a = utils::rev([
            0b0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_u64,
            0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_u64,
            0b0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_u64,
            0b0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0111_u64,
        ]);
        let b = utils::rev([
            0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_u64,
            0b0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_u64,
        ]);

        let mut q = [0u64; 2];
        let mut r = [0u64; 2];

        divrem_bigint(&mut q, Some(&mut r), &mut { a }, &b);
        assert!(all(q == utils::rev([1, 0]), r == utils::rev([0, 7])));
    }

    #[test]
    fn test_div_bigint_1() {
        let a = utils::rev([
            0,
            0b1011101010000111011110010011110100000111110100101000001110111010,
            0b1100010001100101001000000011111111000110110111000011101100100011,
            0b0010011001001000000010000100111100111110101010011101001010110011,
            0b1101101011111000100111001000000101000100110011101110001010011100,
        ]);
        let b = utils::rev([
            0b1010011010101111010001111011000010010000010011100111000110001101,
            0b0111010110010111001111010001100111100100000111101011100000111011,
        ]);

        let mut q = [0; 3];
        let mut r = [0; 2];

        divrem_bigint(&mut q, Some(&mut r), &mut { a }, &b);
        assert!(all(
            q == utils::rev([
                0b1,
                0b0001111001111010010101010110101010011011001111100111011100110001,
                0b0100001010101110101000101110000011100010100111101001010011110001,
            ]),
            r == utils::rev([
                0b1001000110111111111010001101101011111101110101011000101110010101,
                0b1010110101010110101011000000010111110110111101110101011100010001
            ])
        ));
    }

    #[test]
    fn test_div_bigint_3() {
        let a = utils::rev([
            0,
            0b0010011001001000000010000100111100111110101010011101001010110011,
            0b1101101011111000100111001000000101000100110011101110001010011100,
        ]);
        let b = utils::rev([0b1010011010101111010001111011000010010000010011100111000110001101]);

        let mut q = [0u64; 2];
        let mut r = [0u64; 1];

        let a_wide = a[0] as u128 | ((a[1] as u128) << 64);
        let b_wide = b[0] as u128;
        let q_wide = a_wide / b_wide;
        let r_wide = a_wide % b_wide;

        divrem_bigint(&mut q, Some(&mut r), &mut { a }, &b);
        assert!(all(
            q == utils::rev([(q_wide >> 64) as u64, q_wide as u64]),
            r == utils::rev([r_wide as u64])
        ));
    }
}
