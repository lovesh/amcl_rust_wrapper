use crate::constants::GroupG2_SIZE;
use crate::types::GroupG2;
use crate::utils::hash_msg;
use crate::errors::ValueError;
use std::ops::{Add, AddAssign, Sub, Mul, Neg};
use crate::field_elem::{FieldElement};
use crate::group_elem::{GroupElement};

use std::fmt;
use std::slice::Iter;

#[derive(Copy, Clone, Debug)]
pub struct G2 {
    value: GroupG2
}

impl fmt::Display for G2 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let c = self.value.clone();
        write!(f, "{}", c.tostring())
    }
}

impl GroupElement for G2 {
    fn new() -> Self {
        Self {
            value: GroupG2::new()
        }
    }

    fn identity() -> Self {
        let mut v = GroupG2::new();
        v.inf();
        Self {
            value: v
        }
    }

    fn generator() -> Self {
        GroupG2::generator().into()
    }

    fn random() -> Self {
        let n = FieldElement::random();
        Self::generator().scalar_mul_const_time(&n)
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
        self.to_ecp().tostring()
    }
}

/// Represents an element of the sub-group of the elliptic curve over the prime extension field
impl G2 {
    /// Return underlying elliptic curve point, ECP2
    pub fn to_ecp(&self) -> GroupG2 {
        self.value.clone()
    }
}

impl_group_elem_conversions!(G2, GroupG2, GroupG2_SIZE);

impl_group_elem_ops!(G2);

impl_scalar_mul_ops!(G2);

#[cfg(test)]
mod test {
    use super::*;
    use std::time::{Duration, Instant};

    #[test]
    fn test_to_and_from_bytes() {
        for _ in 0..100 {
            let x = G2::random();
            let mut bytes: [u8; GroupG2_SIZE] = [0; GroupG2_SIZE];
            bytes.copy_from_slice(x.to_bytes().as_slice());
            let y = G2::from(&bytes);
            assert_eq!(x, y)
        }
    }

    #[test]
    fn test_negating_group_elems() {
        let b = G2::random();
        let neg_b = -b;
        assert_ne!(b, neg_b);
        let neg_neg_b = -neg_b;
        assert_eq!(b, neg_neg_b);
        assert_eq!(b + neg_b, G2::identity());
    }

    #[test]
    fn test_scalar_mult_operators() {
        for _ in 0..10 {
            let g = G2::random();
            let f = FieldElement::random();
            let m = g.scalar_mul_const_time(&f);
            // Operands can be in any order
            assert_eq!(m, g * f);
            assert_eq!(m, f * g);
        }
    }

    #[test]
    fn test_group_elem_addition() {
        let a = G2::random();
        let b = G2::random();
        let c = G2::random();

        let sum = a + b + c;

        let mut expected_sum = G2::new();
        expected_sum = expected_sum.plus(&a);
        expected_sum = expected_sum.plus(&b);
        expected_sum = expected_sum.plus(&c);
        assert_eq!(sum, expected_sum);
    }

    #[test]
    fn timing_group_elem_addition() {
        let count = 100;
        let points: Vec<_> = (0..100).map(|_| G2::random()).collect();
        let mut R = G2::random();
        let start = Instant::now();
        for i in 0..count {
            R = R + points[i];
        }
        println!("Addition time for {} G2 elems = {:?}", count, start.elapsed());
    }
}