use crate::constants::{CURVE_ORDER, GROUP_G2_SIZE, FIELD_ORDER_ELEMENT_SIZE};
use crate::errors::{SerzDeserzError, ValueError};
use crate::curve_order_elem::{CurveOrderElement, CurveOrderElementVector};
use crate::group_elem::{GroupElement, GroupElementVector};
use crate::types::{GroupG2, FP2, BigNum};
use crate::utils::hash_msg;
use std::iter;
use std::ops::{Add, AddAssign, Index, IndexMut, Mul, Neg, Sub, SubAssign};

use std::fmt;
use std::hash::{Hash, Hasher};
use std::slice::Iter;

use crate::group_elem_g1::parse_hex_as_fp;
use rayon::prelude::*;
use serde::{Serialize, Serializer, Deserialize, Deserializer};
use serde::de::{Error as DError, Visitor};
use std::str::SplitWhitespace;
use zeroize::Zeroize;

#[derive(Clone, Debug)]
pub struct G2 {
    value: GroupG2,
}

impl GroupElement for G2 {
    fn new() -> Self {
        Self {
            value: GroupG2::new(),
        }
    }

    fn identity() -> Self {
        let mut v = GroupG2::new();
        v.inf();
        Self { value: v }
    }

    /// This is an arbitrary choice. Any group element can be a generator
    fn generator() -> Self {
        GroupG2::generator().into()
    }

    fn is_identity(&self) -> bool {
        self.value.is_infinity()
    }

    fn set_to_identity(&mut self) {
        self.value.inf()
    }

    fn from_msg_hash(msg: &[u8]) -> Self {
        GroupG2::mapit(&hash_msg(msg)).into()
    }

    /// TODO: call the appropriate function once implemented in `hash2curve` crate
    fn hash_to_curve(_msg: &[u8], _dst: &hash2curve::DomainSeparationTag) -> Self {
        unimplemented!();
    }

    fn to_vec(&self) -> Vec<u8> {
        let mut bytes: [u8; GROUP_G2_SIZE] = [0; GROUP_G2_SIZE];
        self.write_to_slice_unchecked(&mut bytes);
        bytes.to_vec()
    }

    fn from_slice(bytes: &[u8]) -> Result<Self, SerzDeserzError> {
        if bytes.len() != GROUP_G2_SIZE {
            return Err(SerzDeserzError::G2BytesIncorrectSize(
                bytes.len(),
                GROUP_G2_SIZE,
            ));
        }
        Ok(GroupG2::frombytes(bytes).into())
    }

    fn write_to_slice(&self, target: &mut [u8]) -> Result<(), SerzDeserzError> {
        if target.len() != GROUP_G2_SIZE {
            return Err(SerzDeserzError::G2BytesIncorrectSize(
                target.len(),
                GROUP_G2_SIZE,
            ));
        }
        self.write_to_slice_unchecked(target);
        Ok(())
    }

    fn write_to_slice_unchecked(&self, target: &mut [u8]) {
        let mut temp = GroupG2::new();
        temp.copy(&self.value);
        temp.tobytes(target);
    }

    fn add_assign_(&mut self, b: &Self) {
        self.value.add(&b.value);
    }

    fn sub_assign_(&mut self, b: &Self) {
        self.value.sub(&b.value);
    }

    fn plus(&self, b: &Self) -> Self {
        let mut sum = self.value.clone();
        sum.add(&b.value);
        sum.into()
    }

    fn minus(&self, b: &Self) -> Self {
        let mut diff = self.value.clone();
        diff.sub(&b.value);
        diff.into()
    }

    fn scalar_mul_const_time(&self, a: &CurveOrderElement) -> Self {
        self.value.mul(&a.to_bignum()).into()
    }

    fn double(&self) -> Self {
        let mut d = self.value.clone();
        d.dbl();
        d.into()
    }

    fn double_mut(&mut self) {
        self.value.dbl();
    }

    fn to_hex(&self) -> String {
        self.value.to_hex()
    }

    fn from_hex(s: String) -> Result<Self, SerzDeserzError> {
        let mut iter = s.split_whitespace();
        let x = parse_hex_as_fp2(&mut iter)?;
        let y = parse_hex_as_fp2(&mut iter)?;
        let z = parse_hex_as_fp2(&mut iter)?;
        let mut value = GroupG2::new();
        value.setpx(x);
        value.setpy(y);
        value.setpz(z);
        Ok(Self { value })
    }

    fn negation(&self) -> Self {
        let mut n = self.to_ecp();
        n.neg();
        n.into()
    }

    fn is_extension() -> bool {
        return true;
    }

    fn has_correct_order(&self) -> bool {
        return self.value.mul(&CURVE_ORDER).is_infinity();
    }
}

impl G2 {
    pub fn to_bytes(&self) -> [u8; 4 * FIELD_ORDER_ELEMENT_SIZE] {
        let mut bytes = [0u8; 4 * FIELD_ORDER_ELEMENT_SIZE];
        self.value.tobytes(&mut bytes[..]);
        bytes
    }

    pub fn to_compressed_bytes(&self) -> [u8; 2 * FIELD_ORDER_ELEMENT_SIZE] {
        let mut bytes = [0u8; 2 * FIELD_ORDER_ELEMENT_SIZE];
        let mut temp = GroupG2::new();
        temp.copy(&self.value);
        temp.affine();

        temp.x.geta().tobytes(&mut bytes[..FIELD_ORDER_ELEMENT_SIZE]);
        temp.x.getb().tobytes(&mut bytes[FIELD_ORDER_ELEMENT_SIZE..]);

        let a = temp.y.geta().parity() as u8;
        let b = temp.y.getb().parity() as u8;

        let parity = a << 1 | b;
        bytes[0] |= parity << 6;
        bytes
    }
}

impl From<[u8; 2*FIELD_ORDER_ELEMENT_SIZE]> for G2 {
    fn from(data: [u8; 2*FIELD_ORDER_ELEMENT_SIZE]) -> Self {
        Self::from(&data)
    }
}

impl From<&[u8; 2*FIELD_ORDER_ELEMENT_SIZE]> for G2 {
    fn from(data: &[u8; 2*FIELD_ORDER_ELEMENT_SIZE]) -> Self {
        let mut temp = data.clone();
        let parity = (temp[0] >> 6) & 3u8;
        let pa = if parity & 2u8 == 2 { 1 } else { 0 };
        let pb = if parity & 1u8 == 1 { 1 } else { 0 };
        temp[0] &= 0x3F;
        let mut a = BigNum::frombytes(&temp[..FIELD_ORDER_ELEMENT_SIZE]);
        let mut b = BigNum::frombytes(&temp[FIELD_ORDER_ELEMENT_SIZE..]);
        a.norm();
        b.norm();
        let mut x = FP2::new_bigs(&a, &b);
        x.norm();

        let mut y = GroupG2::rhs(&x);
        if y.sqrt() {
            if y.geta().parity() != pa {
                y.a.neg();
            }
            if y.getb().parity() != pb {
                y.b.neg();
            }
            y.reduce();
        }

        let g2 = GroupG2::new_fp2s(&x, &y);
        Self { value: g2 }
    }
}

/// Parse given hex string as FP2
pub fn parse_hex_as_fp2(iter: &mut SplitWhitespace) -> Result<FP2, SerzDeserzError> {
    // Logic almost copied from AMCL but with error handling and constant time execution.
    // Constant time is important as hex is used during serialization and deserialization.
    // A seemingly effortless solution is to filter string for errors and pad with 0s before
    // passing to AMCL but that would be expensive as the string is scanned twice
    let a = parse_hex_as_fp(iter)?;
    let b = parse_hex_as_fp(iter)?;
    let mut fp2 = FP2::new();
    fp2.seta(a);
    fp2.setb(b);
    Ok(fp2)
}

impl_group_elem_traits!(G2, GroupG2);

impl_group_elem_conversions!(G2, GroupG2, GROUP_G2_SIZE);

impl_group_elem_ops!(G2);

impl_scalar_mul_ops!(G2);

impl_group_element_lookup_table!(G2, G2LookupTable);

// Represents an element of the sub-group of the elliptic curve over prime the extension field
impl_optmz_scalar_mul_ops!(G2, GroupG2, G2LookupTable);

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct G2Vector {
    elems: Vec<G2>,
}

impl_group_elem_vec_ops!(G2, G2Vector);

impl_group_elem_vec_product_ops!(G2, G2Vector, G2LookupTable);

impl_group_elem_vec_conversions!(G2, G2Vector);

impl G2 {
    /// Computes sum of 2 scalar multiplications.
    /// Faster than doing the scalar multiplications individually and then adding them. Uses lookup table
    /// returns self*a + h*b
    pub fn binary_scalar_mul(&self, h: &Self, a: &CurveOrderElement, b: &CurveOrderElement) -> Self {
        // TODO: Replace with faster
        let group_elems = iter::once(self).chain(iter::once(h));
        let field_elems = iter::once(a).chain(iter::once(b));
        G2Vector::multi_scalar_mul_const_time_without_precomputation(group_elems, field_elems)
            .unwrap()
    }
}

#[cfg(test)]
mod test {
    use super::G2;
    use crate::group_elem::GroupElement;
    use crate::curve_order_elem::CurveOrderElement;

    #[test]
    fn test_compression() {
        let g2 = G2::generator();
        let bytes = g2.to_compressed_bytes();
        assert_eq!(G2::from(bytes), g2);

        for _ in 0..30 {
            let sk = CurveOrderElement::random();
            let pk = &g2 * &sk;

            let bytes = pk.to_compressed_bytes();
            let t = G2::from(bytes);

            assert_eq!(t, pk);
        }
    }

    #[test]
    fn test_parse_hex_for_fp2() {
        // TODO:
    }

    #[test]
    fn test_parse_bad_hex_for_fp2() {
        // TODO:
    }
}
