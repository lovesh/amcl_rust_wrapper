use rand::{CryptoRng, RngCore};

use crate::errors::{SerzDeserzError, ValueError};
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

pub trait GroupElement: Clone + Sized {
    fn new() -> Self;

    /// Return the identity element
    fn identity() -> Self;

    /// Return the group's generator
    fn generator() -> Self;

    /// Return a random group element
    fn random() -> Self {
        let n = FieldElement::random();
        Self::generator().scalar_mul_const_time(&n)
    }

    /// Return a random group element using the given random number generator
    fn random_using_rng<R: RngCore + CryptoRng>(rng: &mut R) -> Self {
        let n = FieldElement::random_using_rng(rng);
        Self::generator().scalar_mul_const_time(&n)
    }

    /// Check if the the point is the identity element of the group
    fn is_identity(&self) -> bool;

    /// Set the point to the identity element of the group
    fn set_to_identity(&mut self);

    /// Hash an arbitrary sized message using SHAKE and return output as group element
    fn from_msg_hash(msg: &[u8]) -> Self;

    fn to_bytes(&self) -> Vec<u8>;

    fn from_bytes(bytes: &[u8]) -> Result<Self, SerzDeserzError>;

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

    /// Returns negation of this element
    fn negation(&self) -> Self;

    fn is_extension() -> bool;

    /// Checks if the element has correct order by checking if self *  group order (curve order) == Identity element (point at infinity).
    /// Uses constant time scalar multiplication.
    /// Question: But since we always know the multiplicand (group order) is there a faster way?
    fn has_correct_order(&self) -> bool;

    // TODO: Implement has_correct_order for variable time as well. Need to implement variable time scalar multiplication for group G2.
}

#[macro_export]
macro_rules! impl_group_elem_conversions {
    ( $group_element:ident, $group:ident, $group_size:ident ) => {
        impl From<$group> for $group_element {
            fn from(x: $group) -> Self {
                Self { value: x }
            }
        }

        impl From<&$group> for $group_element {
            fn from(x: &$group) -> Self {
                Self { value: x.clone() }
            }
        }

        impl From<&[u8; $group_size]> for $group_element {
            fn from(x: &[u8; $group_size]) -> Self {
                Self {
                    value: $group::frombytes(x),
                }
            }
        }

        impl Hash for $group_element {
            fn hash<H: Hasher>(&self, state: &mut H) {
                state.write(&self.to_bytes())
            }
        }
    };
}

#[macro_export]
macro_rules! impl_group_elem_traits {
    ( $group_element:ident ) => {
        impl Default for $group_element {
            fn default() -> Self { Self::new() }
        }

        impl fmt::Display for $group_element {
            fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
                let c = self.value.clone();
                write!(f, "{}", c.tostring())
            }
        }
    };
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

        impl Eq for $group_element {}

        impl Add for $group_element {
            type Output = Self;

            fn add(self, other: Self) -> Self {
                self.plus(&other)
            }
        }

        impl Add<$group_element> for &$group_element {
            type Output = $group_element;

            fn add(self, other: $group_element) -> $group_element {
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
    };
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
    fn plus(&self, b: &Self) -> Result<Self, ValueError>;

    /// Subtract 2 vectors
    fn minus(&self, b: &Self) -> Result<Self, ValueError>;

    fn iter(&self) -> Iter<T>;
}

#[macro_export]
macro_rules! impl_group_elem_vec_ops {
    ( $group_element:ident, $group_element_vec:ident ) => {
        impl From<Vec<$group_element>> for $group_element_vec {
            fn from(x: Vec<$group_element>) -> Self {
                Self { elems: x }
            }
        }

        impl From<&[$group_element]> for $group_element_vec {
            fn from(x: &[$group_element]) -> Self {
                Self { elems: x.to_vec() }
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
                    return false;
                }
                for i in 0..self.len() {
                    if self[i] != other[i] {
                        return false;
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
    };
}

#[cfg(test)]
mod test {
    use super::*;
    use std::collections::{HashMap, HashSet};
    use std::time::{Duration, Instant};
    use crate::group_elem_g1::G1;
    use crate::group_elem_g2::G2;
    use crate::constants::{GroupG1_SIZE, GroupG2_SIZE};

    #[test]
    fn test_to_and_from_bytes() {
        macro_rules! to_and_fro_bytes {
            ( $group:ident, $group_size:ident ) => {
                let x = $group::random();
                let mut bytes: [u8; $group_size] = [0; $group_size];
                bytes.copy_from_slice(x.to_bytes().as_slice());
                let y = $group::from(&bytes);
                assert_eq!(x, y);

                let bytes1 = x.to_bytes();
                assert_eq!(x, $group::from_bytes(bytes1.as_slice()).unwrap());

                // Increase length of byte vector by adding a byte. Choice of byte is arbitrary
                let mut bytes2 = bytes1.clone();
                bytes2.push(0);
                assert!($group::from_bytes(bytes2.as_slice()).is_err());

                // Decrease length of byte vector
                assert!($group::from_bytes(&bytes1[0..$group_size - 1]).is_err());
            };
        }

        let count = 100;
        for _ in 0..count {
            to_and_fro_bytes!(G1, GroupG1_SIZE);
        }

        for _ in 0..count {
            to_and_fro_bytes!(G2, GroupG2_SIZE);
        }
    }

    #[test]
    fn test_hashing() {
        // If the element can be added to HashSet or HashMap, it must be hashable.
        macro_rules! hashing {
            ( $group:ident ) => {
                {
                    let mut set = HashSet::new();
                    let mut map = HashMap::new();
                    set.insert($group::random());
                    map.insert($group::random(), $group::random());
                }
            };
        }

        hashing!(G1);
        hashing!(G2);
    }

    #[test]
    fn test_negating_group_elems() {
        macro_rules! negating {
            ( $group:ident ) => {
                {
                    let b = $group::random();
                    let neg_b = -b;
                    assert_ne!(b, neg_b);
                    let neg_neg_b = -neg_b;
                    assert_eq!(b, neg_neg_b);
                    assert_eq!(b + neg_b, $group::identity());
                }
            };
        }
        negating!(G1);
        negating!(G2);
    }

    #[test]
    fn test_scalar_mult_operators() {
        macro_rules! scalar_mult {
            ( $group:ident ) => {
                {
                    let g = $group::random();
                    let f = FieldElement::random();
                    let m = g.scalar_mul_const_time(&f);
                    // Operands can be in any order
                    assert_eq!(m, &g * &f);
                    assert_eq!(m, &f * &g);
                }
            };
        }

        for _ in 0..10 {
            scalar_mult!(G1)
        }
        for _ in 0..10 {
            scalar_mult!(G2)
        }
    }

    #[test]
    fn test_group_elem_addition() {
        macro_rules! addition {
            ( $group:ident ) => {
                {
                    let a = G1::random();
                    let b = G1::random();
                    let c = G1::random();

                    let sum = a + b + c;

                    let mut expected_sum = G1::new();
                    expected_sum = expected_sum.plus(&a);
                    expected_sum = expected_sum.plus(&b);
                    expected_sum = expected_sum.plus(&c);
                    assert_eq!(sum, expected_sum);
                }
            };
        }
        addition!(G1);
        addition!(G2);
    }

    #[test]
    fn test_negation() {
        macro_rules! neg {
            ( $group:ident ) => {
                {
                    let a = G1::random();
                    let b = a.negation();
                    assert!((a + b).is_identity())
                }
            };
        }

        for i in 0..10 {
            neg!(G1);
            neg!(G2);
        }
    }

    #[test]
    fn timing_correct_order_check() {
        let count = 10;
        macro_rules! order_check {
            ( $group:ident ) => {
                {
                    let start = Instant::now();
                    for _ in 0..count {
                        let a = $group::random();
                        assert!(a.has_correct_order())
                    }
                    println!(
                        "For {} elements, time to check correct order is {:?}",
                        count,
                        start.elapsed()
                    )
                }
            };
        }
        order_check!(G1);
        order_check!(G2);
    }

    #[test]
    fn timing_group_elem_addition_and_scalar_multiplication() {
        let count = 100;
        macro_rules! add_mul {
            ( $group:ident ) => {
                {
                    let points: Vec<_> = (0..100).map(|_| $group::random()).collect();
                    let mut R = $group::random();
                    let mut start = Instant::now();
                    for i in 0..count {
                        R = R + points[i];
                    }
                    println!(
                        "Addition time for {} elems = {:?}",
                        count,
                        start.elapsed()
                    );

                    let fs: Vec<_> = (0..100).map(|_| FieldElement::random()).collect();
                    start = Instant::now();
                    for i in 0..count {
                        points[i] * &fs[i];
                    }
                    println!(
                        "Scalar multiplication time for {} elems = {:?}",
                        count,
                        start.elapsed()
                    );
                }
            };
        }

        add_mul!(G1);
        add_mul!(G2);
    }
}
