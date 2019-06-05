use super::ECCurve::rom;
use super::ECCurve::big::{NLEN as curve_NLEN, MODBYTES as curve_MODBYTES};
use super::types::{BigNum, DoubleBigNum, GroupG1, GroupG2};

pub const MODBYTES: usize = curve_MODBYTES;
pub const NLEN: usize = curve_NLEN;

// TODO: Extra 1 byte not needed.
// Byte size of element in group G1, 1 extra byte for compression flag
pub const GroupG1_SIZE: usize = (2 * MODBYTES + 1) as usize;
// Byte size of element in group G2
pub const GroupG2_SIZE: usize = (4 * MODBYTES) as usize;

lazy_static! {
    pub static ref GeneratorG1: GroupG1 = GroupG1::generator();
    pub static ref GeneratorG2: GroupG2 = GroupG2::generator();
    pub static ref CurveOrder: BigNum = BigNum::new_ints(&rom::CURVE_ORDER);
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
