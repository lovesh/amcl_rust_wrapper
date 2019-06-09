use crate::errors::{ValueError, SerzDeserzError};
use crate::field_elem::FieldElement;
use std::slice::Iter;

#[macro_export]
macro_rules! add_group_elems {
    ( $( $elem:expr ),* ) => {
        {
            let mut sum = GroupElement::new();
            $(
                sum += $elem;
            )*
            sum
        }
    };
}

pub trait GroupElement: Sized {
    fn new() -> Self;

    /// Return the identity element
    fn identity() -> Self;

    /// Return the group's generator
    fn generator() -> Self;

    fn random() -> Self {
        let n = FieldElement::random();
        Self::generator().scalar_mul_const_time(&n)
    }

    /// Check if the the point is the identity element of the group
    fn is_identity(&self) -> bool;

    /// Set the point to the identity element of the group
    fn set_to_identity(&mut self);

    /// Hash an arbitrary sized message using SHAKE and return output as group element
    fn from_msg_hash(msg: &[u8]) -> Self;

    fn to_bytes(&self) -> Vec<u8>;

    fn from_bytes(bytes: &[u8])  -> Result<Self, SerzDeserzError>;

    /// Add a group element to itself. `self = self + b`
    fn add_assign_(&mut self, b: &Self);

    /// Subtract a group element from itself. `self = self - b`
    fn sub_assign_(&mut self, b: &Self);

    /// Return sum of a group element and itself. `self + b`
    fn plus(&self, b: &Self) -> Self;

    /// Return difference of a group element and itself. `self - b`
    fn minus(&self, b: &Self) -> Self;

    /// Multiply point on the curve (element of group G1) with a scalar. Constant time operation.
    /// self * field_element_a.
    fn scalar_mul_const_time(&self, a: &FieldElement) -> Self;

    /// Return the double of the group element
    fn double(&self) -> Self;

    fn double_mut(&mut self);

    /// Returns hex string
    fn to_hex(&self) -> String;
}

#[macro_export]
macro_rules! impl_group_elem_conversions {
    ( $group_element:ident, $group:ident, $group_size:ident ) => {
        impl From<$group> for $group_element {
            fn from(x: $group) -> Self {
                Self {
                    value: x
                }
            }
        }

        impl From<&$group> for $group_element {
            fn from(x: &$group) -> Self {
                Self {
                value: x.clone()
                }
            }
        }

        impl From<&[u8; $group_size]> for $group_element {
            fn from(x: &[u8; $group_size]) -> Self {
                Self {
                    value: $group::frombytes(x)
                }
            }
        }
    }
}

#[macro_export]
macro_rules! impl_group_elem_ops {
    ( $group_element:ident ) => {
        impl PartialEq for $group_element {
            fn eq(&self, other: &$group_element) -> bool {
                let l = self.clone();
                let mut r = other.clone();
                l.value.equals(&mut r.value)
            }
        }

        impl Add for $group_element {
            type Output = Self;

            fn add(self, other: Self) -> Self {
                self.plus(&other)
            }
        }

        impl Add<$group_element> for &$group_element {
            type Output = $group_element;

            fn add(self, other: $group_element) ->$group_element {
                self.plus(&other)
            }
        }

        impl<'a> Add<&'a $group_element> for $group_element {
            type Output = Self;
            fn add(self, other: &'a $group_element) -> Self {
                self.plus(other)
            }
        }

        impl AddAssign for $group_element {
            fn add_assign(&mut self, other: Self) {
                self.add_assign_(&other)
            }
        }

        impl Sub for $group_element {
            type Output = Self;

            fn sub(self, other: Self) -> Self {
                self.minus(&other)
            }
        }

        impl Mul<FieldElement> for $group_element {
            type Output = Self;

            fn mul(self, other: FieldElement) -> Self {
                self.scalar_mul_const_time(&other)
            }
        }

        impl Mul<&FieldElement> for $group_element {
            type Output = Self;

            fn mul(self, other: &FieldElement) -> Self {
                self.scalar_mul_const_time(other)
            }
        }

        impl Mul<FieldElement> for &$group_element {
            type Output = $group_element;

            fn mul(self, other: FieldElement) -> $group_element {
                self.scalar_mul_const_time(&other)
            }
        }

        impl Mul<&FieldElement> for &$group_element {
            type Output = $group_element;

            fn mul(self, other: &FieldElement) -> $group_element {
                self.scalar_mul_const_time(other)
            }
        }

        impl Neg for $group_element {
            type Output = Self;

            fn neg(self) -> Self::Output {
                let mut t = self.to_ecp();
                t.neg();
                t.into()
            }
        }

        impl Neg for &$group_element {
            type Output = $group_element;

            fn neg(self) -> Self::Output {
                let mut t = self.to_ecp();
                t.neg();
                t.into()
            }
        }
    }
}

macro_rules! impl_scalar_mul_ops {
    ( $group_element:ident ) => {

        impl Mul<$group_element> for FieldElement {
            type Output = $group_element;

            fn mul(self, other: $group_element) -> $group_element {
                other.scalar_mul_const_time(&self)
            }
        }

        impl Mul<&$group_element> for FieldElement {
            type Output = $group_element;

            fn mul(self, other: &$group_element) -> $group_element {
                other.scalar_mul_const_time(&self)
            }
        }

        impl Mul<$group_element> for &FieldElement {
            type Output = $group_element;

            fn mul(self, other: $group_element) -> $group_element {
                other.scalar_mul_const_time(self)
            }
        }

        impl Mul<&$group_element> for &FieldElement {
            type Output = $group_element;

            fn mul(self, other: &$group_element) -> $group_element {
                other.scalar_mul_const_time(self)
            }
        }
    };
}

pub trait GroupElementVector<T>: Sized {
    fn new(size: usize) -> Self;

    fn with_capacity(capacity: usize) -> Self;

    fn as_slice(&self) -> &[T];

    fn len(&self) -> usize;

    fn push(&mut self, value: T);

    fn append(&mut self, other: &mut Self);

    /// Compute sum of all elements of a vector
    fn sum(&self) -> T;

    /// Multiply each field element of the vector with another given field
    /// element `n` (scale the vector)
    fn scale(&mut self, n: &FieldElement);

    fn scaled_by(&self, n: &FieldElement) -> Self;

    /// Add 2 vectors
    fn plus(&self, b: &Self) ->  Result<Self, ValueError>;

    /// Subtract 2 vectors
    fn minus(&self, b: &Self) ->  Result<Self, ValueError>;

    fn iter(&self) -> Iter<T>;
}

#[macro_export]
macro_rules! impl_group_elem_vec_ops {
    ( $group_element:ident, $group_element_vec:ident ) => {
        impl From<Vec<$group_element>> for $group_element_vec {
            fn from(x: Vec<$group_element>) -> Self {
                Self {
                    elems: x
                }
            }
        }

        impl From<&[$group_element]> for $group_element_vec {
            fn from(x: &[$group_element]) -> Self {
                Self {
                    elems: x.to_vec()
                }
            }
        }

        impl Index<usize> for $group_element_vec {
            type Output = $group_element;

            fn index(&self, idx: usize) -> &$group_element {
                &self.elems[idx]
            }
        }

        impl IndexMut<usize> for $group_element_vec {

            fn index_mut(&mut self, idx: usize) -> &mut $group_element {
                &mut self.elems[idx]
            }
        }

        impl PartialEq for $group_element_vec {
            fn eq(&self, other: &Self) -> bool {
                if self.len() != other.len() {
                    return false
                }
                for i in 0..self.len() {
                    if self[i] != other[i] {
                        return false
                    }
                }
                true
            }
        }

        impl IntoIterator for $group_element_vec {
            type Item = $group_element;
            type IntoIter = ::std::vec::IntoIter<$group_element>;

            fn into_iter(self) -> Self::IntoIter {
                self.elems.into_iter()
            }
        }
    }
}