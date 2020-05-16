use crate::amcl::hmac;
use crate::constants::{CurveOrder, GroupG2_SIZE, G2_COMP_BYTE_SIZE, HASH_TYPE};
use crate::errors::{SerzDeserzError, ValueError};
use crate::field_elem::{FieldElement, FieldElementVector};
use crate::group_elem::{GroupElement, GroupElementVector};
use crate::types::{GroupG2, FP, FP2};
use crate::utils::{hash_msg, hash_to_field};
use std::iter;
use std::ops::{Add, AddAssign, Index, IndexMut, Mul, Neg, Sub, SubAssign};

use core::fmt;
use std::hash::{Hash, Hasher};
use std::slice::Iter;

use crate::rayon::iter::IntoParallelRefMutIterator;
use rayon::prelude::*;
use serde::de::{Deserialize, Deserializer, Error as DError, Visitor};
use serde::ser::{Error as SError, Serialize, Serializer};
use zeroize::Zeroize;

/// Don't derive Copy trait as it can hold secret data and should not be accidentally copied
#[derive(Clone)]
pub struct G2 {
    value: GroupG2,
}

impl fmt::Debug for G2 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "ECP2: [ {} ]", self.value.tostring())
    }
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

    impl_group_elem_byte_conversion_methods!(
        GroupG2,
        GroupG2_SIZE,
        G2_COMP_BYTE_SIZE,
        SerzDeserzError::G2BytesIncorrectSize
    );

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

    /// Returns the string `infinity` if the element corresponds to a point at infinity
    /// Returns `(x,y)` where both `x` and `y` are hex representations of FP2
    fn to_hex(&self) -> String {
        self.value.tostring()
    }

    fn from_hex(mut string: String) -> Result<Self, SerzDeserzError> {
        if &string == "infinity" {
            return Ok(Self::new());
        }

        // Need string as "(x,y)"
        unbound_bounded_string!(string, '(', ')', SerzDeserzError::CannotParseG2);

        let (x, y) = split_string_to_2_tuple!(string, SerzDeserzError::CannotParseG2);

        let x_fp2 = parse_hex_as_FP2(x)?;
        let y_fp2 = parse_hex_as_FP2(y)?;

        Ok(Self {
            value: GroupG2::new_fp2s(&x_fp2, &y_fp2),
        })
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
pub fn parse_hex_as_FP2(mut string: String) -> Result<FP2, SerzDeserzError> {
    // Need string as "[a,b]"
    unbound_bounded_string!(string, '[', ']', SerzDeserzError::CannotParseFP2);

    let (a, b) = split_string_to_2_tuple!(string, SerzDeserzError::CannotParseFP2);

    let a_big = FieldElement::parse_hex_as_bignum(a)?;
    let b_big = FieldElement::parse_hex_as_bignum(b)?;
    Ok(FP2::new_bigs(&a_big, &b_big))
}

impl_group_elem_traits!(G2, GroupG2);

impl_group_elem_serz!(G2, GroupG2, "G2");

impl_group_elem_conversions!(G2, GroupG2, GroupG2_SIZE, G2_COMP_BYTE_SIZE);

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

impl G2 {
    /// Computes sum of 2 scalar multiplications.
    /// Faster than doing the scalar multiplications individually and then adding them. Uses lookup table
    /// returns self*a + h*b
    pub fn binary_scalar_mul(&self, h: &Self, a: &FieldElement, b: &FieldElement) -> Self {
        // TODO: Replace with faster
        let group_elems = iter::once(self).chain(iter::once(h));
        let field_elems = iter::once(a).chain(iter::once(b));
        G2Vector::multi_scalar_mul_const_time_without_precomputation(group_elems, field_elems)
            .unwrap()
    }

    /// Hashes a byte slice to a group element according to the hash to curve point IETF standard
    /// https://datatracker.ietf.org/doc/draft-irtf-cfrg-hash-to-curve/?include_text=1
    /// `domain_separation_tag` should be unique between protocols as well as curves, eg. protocol A and
    /// protocol B should use different `domain_separation_tag` while hashing to the same curve and
    /// protocol A should use different `domain_separation_tag` while hashing to different curves.
    /// Look at section 3.1 of the standard for more details
    pub fn hash_to_curve(dst: &[u8], msg: &[u8]) -> G2 {
        // Get 4 field elements as FP
        let mut u: [FP; 4] = [FP::new(), FP::new(), FP::new(), FP::new()];
        hash_to_field(hmac::MC_SHA2, HASH_TYPE, dst, msg, &mut u, 4);

        // Create extension field elements (FP^2) as FP2
        let fp2_1 = FP2::new_fps(&u[0], &u[1]);
        let fp2_2 = FP2::new_fps(&u[2], &u[3]);

        // Map each FP2 to a curve point and add the points
        let mut P = GroupG2::map2point(&fp2_1);
        let P1 = GroupG2::map2point(&fp2_2);
        P.add(&P1);
        // clear the cofactor of the addition point
        P.cfp();

        Self { value: P }
    }
}

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
