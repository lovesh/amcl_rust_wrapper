use super::types::{BigNum, DoubleBigNum, GroupG1};

use super::ECCurve::big::{BASEBITS, MODBYTES as curve_MODBYTES, NLEN as curve_NLEN};
use super::ECCurve::rom;

pub const MODBYTES: usize = curve_MODBYTES;
pub const NLEN: usize = curve_NLEN;
pub const BigNumBits: usize = BASEBITS;

// Byte size of element in group G1, 1 extra byte for compression flag
pub const FieldElement_SIZE: usize = MODBYTES;

// Byte size of uncompressed element in group G1, 1 extra byte for compression flag
pub const GroupG1_SIZE: usize = (2 * MODBYTES + 1) as usize;

// Byte size of compressed element in group G1, 1 extra byte for compression flag
pub const GroupG1_COMP_SIZE: usize = (MODBYTES + 1) as usize;

lazy_static! {
    pub static ref GeneratorG1: GroupG1 = GroupG1::generator();
    pub static ref CurveOrder: BigNum = BigNum::new_ints(&rom::CURVE_ORDER);
    pub static ref CurveOrderBitSize: usize = CurveOrder.nbits();
    pub static ref FieldElementZero: BigNum = BigNum::new();
    pub static ref BarrettRedc_k: usize = CurveOrder.nbits();
    pub static ref BarrettRedc_u: BigNum = {
        let k = CurveOrder.nbits();
        let mut u = DoubleBigNum::new();
        u.w[0] = 1;
        // `u.shl(2*k)` crashes, so perform shl(k) twice
        u.shl(k);
        u.shl(k);

        // div returns floored value
        u.div(&CurveOrder)
    };

    pub static ref BarrettRedc_v: BigNum = {
        let k = CurveOrder.nbits();
        let mut v = BigNum::new_int(1isize);
        v.shl(k+1);
        v
    };
}

#[cfg(any(feature = "bls381", feature = "bn254"))]
pub use crate::types_g2::{GeneratorG2, GroupG2_SIZE, GroupGT_SIZE};
