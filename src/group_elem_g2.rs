use crate::constants::{CurveOrder, GroupG2_SIZE};
use crate::errors::{SerzDeserzError, ValueError};
use crate::field_elem::FieldElement;
use crate::group_elem::{GroupElement, GroupElementVector};
use crate::types::GroupG2;
use crate::utils::hash_msg;
use std::ops::{Add, AddAssign, Index, IndexMut, Mul, Neg, Sub};

use std::fmt;
use std::hash::{Hash, Hasher};
use std::slice::Iter;

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
        self.to_ecp().tostring()
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

/// Represents an element of the sub-group of the elliptic curve over the prime extension field
impl G2 {
    /// Return underlying elliptic curve point, ECP2
    pub fn to_ecp(&self) -> GroupG2 {
        self.value.clone()
    }
}

impl_group_elem_traits!(G2);

impl_group_elem_conversions!(G2, GroupG2, GroupG2_SIZE);

impl_group_elem_ops!(G2);

impl_scalar_mul_ops!(G2);

#[derive(Clone, Debug)]
pub struct G2Vector {
    elems: Vec<G2>,
}

impl GroupElementVector<G2> for G2Vector {
    fn new(size: usize) -> Self {
        Self {
            elems: (0..size).map(|_| G2::new()).collect(),
        }
    }

    fn with_capacity(capacity: usize) -> Self {
        Self {
            elems: Vec::<G2>::with_capacity(capacity),
        }
    }

    fn as_slice(&self) -> &[G2] {
        &self.elems
    }

    fn len(&self) -> usize {
        self.elems.len()
    }

    fn push(&mut self, value: G2) {
        self.elems.push(value)
    }

    fn append(&mut self, other: &mut Self) {
        self.elems.append(&mut other.elems)
    }

    fn sum(&self) -> G2 {
        let mut accum = G2::new();
        for i in 0..self.len() {
            accum += self[i];
        }
        accum
    }

    fn scale(&mut self, n: &FieldElement) {
        for i in 0..self.len() {
            self[i] = self[i] * n;
        }
    }

    fn scaled_by(&self, n: &FieldElement) -> Self {
        let mut scaled = Self::with_capacity(self.len());
        for i in 0..self.len() {
            scaled.push(self[i] * n)
        }
        scaled.into()
    }

    fn plus(&self, b: &Self) -> Result<Self, ValueError> {
        check_vector_size_for_equality!(self, b)?;
        let mut sum_vector = Self::with_capacity(self.len());
        for i in 0..self.len() {
            sum_vector.push(self[i] + b.elems[i])
        }
        Ok(sum_vector)
    }

    fn minus(&self, b: &Self) -> Result<Self, ValueError> {
        check_vector_size_for_equality!(self, b)?;
        let mut diff_vector = Self::with_capacity(self.len());
        for i in 0..self.len() {
            diff_vector.push(self[i] - b[i])
        }
        Ok(diff_vector)
    }

    fn iter(&self) -> Iter<G2> {
        self.as_slice().iter()
    }
}

impl_group_elem_vec_ops!(G2, G2Vector);
