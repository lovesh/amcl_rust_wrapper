use super::types::{BigNum, DoubleBigNum, GroupG1};

use super::ECCurve::big::{BASEBITS, MODBYTES as curve_MODBYTES, NLEN as curve_NLEN};
use super::ECCurve::rom;

pub const MODBYTES: usize = curve_MODBYTES;
pub const NLEN: usize = curve_NLEN;
pub const BIG_NUM_BITS: usize = BASEBITS;

// Byte size of element in group G1, 1 extra byte for compression flag
pub const FIELD_ELEMENT_SIZE: usize = MODBYTES;

// Byte size of element in group G1, 1 extra byte for compression flag
pub const GROUP_G1_SIZE: usize = (2 * MODBYTES + 1) as usize;

lazy_static! {
    pub static ref GENERATOR_G1: GroupG1 = GroupG1::generator();
    pub static ref CURVE_ORDER: BigNum = BigNum::new_ints(&rom::CURVE_ORDER);
    pub static ref CURVE_ORDER_BIT_SIZE: usize = CURVE_ORDER.nbits();
    pub static ref FIELD_ELEMENT_ZERO: BigNum = BigNum::new();
    pub static ref BARRETT_REDC_K: usize = CURVE_ORDER.nbits();
    pub static ref BARRETT_REDC_U: BigNum = {
        let k = CURVE_ORDER.nbits();
        let mut u = DoubleBigNum::new();
        u.w[0] = 1;
        // `u.shl(2*k)` crashes, so perform shl(k) twice
        u.shl(k);
        u.shl(k);

        // div returns floored value
        u.div(&CURVE_ORDER)
    };

    pub static ref BARRETT_REDC_V: BigNum = {
        let k = CURVE_ORDER.nbits();
        let mut v = BigNum::new_int(1isize);
        v.shl(k+1);
        v
    };
}

#[cfg(any(feature = "bls381", feature = "bn254"))]
pub use crate::types_g2::{GENERATOR_G2, GROUP_G2_SIZE, GROUP_GT_SIZE};
