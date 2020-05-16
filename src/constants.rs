use super::types::{BigNum, DoubleBigNum};

use super::ECCurve::big::{
    BASEBITS, DNLEN as curve_DNLEN, MODBYTES as curve_MODBYTES, NLEN as curve_NLEN,
};
use super::ECCurve::rom;

pub use super::ECCurve::ecp::{AESKEY, CURVETYPE, EDWARDS, HASH_TYPE, MONTGOMERY, WEIERSTRASS};
pub const MODBYTES: usize = curve_MODBYTES;
pub const NLEN: usize = curve_NLEN;
pub const DNLEN: usize = curve_NLEN;
pub const BigNumBits: usize = BASEBITS;

// Byte size of element in group G1, 1 extra byte for compression flag
pub const FieldElement_SIZE: usize = MODBYTES;

// TODO: For Montgomery curves, the byte size is always MODBYTES. Handle this with features
// Byte size of element in group G1, 1 extra byte for compression flag
pub const GroupG1_SIZE: usize = 2 * MODBYTES + 1;
pub const G1_COMP_BYTE_SIZE: usize = MODBYTES + 1;

fn ceil(a: &usize, b: &usize) -> usize {
    return (a - 1) / b + 1;
}

lazy_static! {
    pub static ref CurveOrder: BigNum = BigNum::new_ints(&rom::CURVE_ORDER);
    pub static ref CurveOrderBitSize: usize = CurveOrder.nbits();
    pub static ref FieldModulus: BigNum = BigNum::new_ints(&rom::MODULUS);
    pub static ref FieldModulusBitSize: usize = FieldModulus.nbits();
    // Used in hashing arbitrary message to curve point
    pub static ref Ell: usize = ceil(&(*FieldModulusBitSize + AESKEY * 8), &8);

    // Constants for Barrett reduction for reducing modulo CurveOrder
    pub static ref BarrettRedc_k: usize = *CurveOrderBitSize;
    pub static ref BarrettRedc_u: BigNum = {
        let k = *BarrettRedc_k;
        let mut u = DoubleBigNum::new();
        u.w[0] = 1;
        // `u.shl(2*k)` crashes, so perform shl(k) twice
        u.shl(k);
        u.shl(k);

        // div returns floored value
        u.div(&CurveOrder)
    };
    pub static ref BarrettRedc_v: BigNum = {
        let k = *BarrettRedc_k;
        let mut v = BigNum::new_int(1isize);
        v.shl(k+1);
        v
    };

    // Constants for Barrett reduction for reducing modulo FieldModulus
    pub static ref BarrettRedc_FM_k: usize = *FieldModulusBitSize;
    pub static ref BarrettRedc_FM_u: BigNum = {
        let k = *BarrettRedc_FM_k;
        let mut u = DoubleBigNum::new();
        u.w[0] = 1;
        // `u.shl(2*k)` crashes, so perform shl(k) twice
        u.shl(k);
        u.shl(k);

        // div returns floored value
        u.div(&FieldModulus)
    };
    pub static ref BarrettRedc_FM_v: BigNum = {
        let k = *BarrettRedc_FM_k;
        let mut v = BigNum::new_int(1isize);
        v.shl(k+1);
        v
    };
}

#[cfg(any(feature = "bls381", feature = "bn254"))]
pub use crate::types_g2::{GroupG2_SIZE, GroupGT_SIZE, G2_COMP_BYTE_SIZE};
