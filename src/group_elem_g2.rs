use clear_on_drop::clear::Clear;

use crate::constants::{CurveOrder, GroupG2_SIZE};
use crate::errors::{SerzDeserzError, ValueError};
use crate::field_elem::{FieldElement, FieldElementVector};
use crate::group_elem::{GroupElement, GroupElementVector};
use crate::types::{GroupG2, FP2};
use crate::utils::hash_msg;
use std::ops::{Add, AddAssign, Index, IndexMut, Mul, Neg, Sub};

use std::fmt;
use std::hash::{Hash, Hasher};
use std::slice::Iter;

use crate::group_elem_g1::parse_hex_as_FP;
use serde::de::{Deserialize, Deserializer, Error as DError, Visitor};
use serde::ser::{Error as SError, Serialize, Serializer};
use std::str::SplitWhitespace;

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

    fn to_bytes(&self) -> Vec<u8> {
        let mut temp = GroupG2::new();
        temp.copy(&self.value);
        let mut bytes: [u8; GroupG2_SIZE] = [0; GroupG2_SIZE];
        temp.tobytes(&mut bytes);
        bytes.to_vec()
    }

    fn from_bytes(bytes: &[u8]) -> Result<Self, SerzDeserzError> {
        if bytes.len() != GroupG2_SIZE {
            return Err(SerzDeserzError::G2BytesIncorrectSize(
                bytes.len(),
                GroupG2_SIZE,
            ));
        }
        Ok(GroupG2::frombytes(bytes).into())
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

    fn scalar_mul_const_time(&self, a: &FieldElement) -> Self {
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
        let x = parse_hex_as_FP2(&mut iter)?;
        let y = parse_hex_as_FP2(&mut iter)?;
        let z = parse_hex_as_FP2(&mut iter)?;
        let mut value = GroupG2::new();
        value.setpx(x);
        value.setpy(y);
        value.setpz(z);
        Ok(G2 { value })
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
        return self.value.mul(&CurveOrder).is_infinity();
    }
}

/// Parse given hex string as FP2
pub fn parse_hex_as_FP2(iter: &mut SplitWhitespace) -> Result<FP2, SerzDeserzError> {
    // Logic almost copied from AMCL but with error handling and constant time execution.
    // Constant time is important as hex is used during serialization and deserialization.
    // A seemingly effortless solution is to filter string for errors and pad with 0s before
    // passing to AMCL but that would be expensive as the string is scanned twice
    let a = parse_hex_as_FP(iter)?;
    let b = parse_hex_as_FP(iter)?;
    let mut fp2 = FP2::new();
    fp2.seta(a);
    fp2.setb(b);
    Ok(fp2)
}

impl_group_elem_traits!(G2);

impl_group_elem_conversions!(G2, GroupG2, GroupG2_SIZE);

impl_group_elem_ops!(G2);

impl_scalar_mul_ops!(G2);

impl_group_element_lookup_table!(G2, G2LookupTable);

/// Represents an element of the sub-group of the elliptic curve over prime the extension field
impl_optmz_scalar_mul_ops!(G2, GroupG2, G2LookupTable);

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct G2Vector {
    elems: Vec<G2>,
}

impl_group_elem_vec_ops!(G2, G2Vector);

impl_group_elem_vec_product_ops!(G2, G2Vector, G2LookupTable);

impl_group_elem_vec_conversions!(G2, G2Vector);

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_parse_hex_for_FP2() {
        // TODO:
    }

    #[test]
    fn test_parse_bad_hex_for_FP2() {
        // TODO:
    }
}
