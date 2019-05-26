use crate::constants::{MODBYTES, GroupG1_SIZE};
use crate::types::{BigNum, GroupG1};
#[macro_use]
use crate::macros;
use crate::utils::hash_msg;
use crate::errors::ValueError;
use std::cmp::Ordering;
use std::ops::{Index, IndexMut, Add, AddAssign, Sub, Mul};
use crate::field_elem::{FieldElement, FieldElementVector};
use std::fmt;
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

#[derive(Copy, Clone, Debug)]
pub struct GroupElement {
    value: GroupG1
}

impl fmt::Display for GroupElement {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let c = self.value.clone();
        write!(f, "{}", c.tostring())
    }
}

/// Represents an element of the group on the elliptic curve
impl GroupElement {
    pub fn new() -> Self {
        Self {
            value: GroupG1::new()
        }
    }

    /// Return the identity element
    pub fn identity() -> Self {
        let mut v = GroupG1::new();
        v.inf();
        Self {
            value: v
        }
    }

    pub fn generator() -> Self {
        GroupG1::generator().into()
    }

    pub fn random(order: Option<&BigNum>) -> Self {
        let n = FieldElement::random(order);
        Self::generator().scalar_mul_const_time(&n)
    }

    pub fn is_identity(&self) -> bool {
        self.value.is_infinity()
    }

    pub fn set_to_identity(&mut self) {
        self.value.inf()
    }

    /// Hash message and return output as group element
    pub fn from_msg_hash(msg: &[u8]) -> GroupElement {
        GroupG1::mapit(&hash_msg(msg)).into()
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut temp = GroupG1::new();
        temp.copy(&self.value);
        let mut bytes: [u8; GroupG1_SIZE] = [0; GroupG1_SIZE];
        temp.tobytes(&mut bytes, false);
        bytes.to_vec()
    }

    pub fn to_ecp(&self) -> GroupG1 {
        self.value.clone()
    }

    /// Add a group element to itself. `self = self + b`
    pub fn add_assign_(&mut self, b: &Self) {
        self.value.add(&b.value);
    }

    /// Subtract a group element from itself. `self = self - b`
    pub fn sub_assign_(&mut self, b: &Self) {
        self.value.sub(&b.value);
    }

    /// Return sum of a group element and itself. `self + b`
    pub fn plus(&self, b: &Self) -> Self {
        let mut sum = self.value.clone();
        sum.add(&b.value);
        sum.into()
    }

    /// Return difference of a group element and itself. `self - b`
    pub fn minus(&self, b: &Self) -> Self {
        let mut diff = self.value.clone();
        diff.sub(&b.value);
        diff.into()
    }

    /// Multiply point on the curve (element of group G1) with a scalar. Constant time operation.
    /// self * field_element_a.
    pub fn scalar_mul_const_time(&self, a: &FieldElement) -> Self {
        self.value.mul(&a.to_bignum()).into()
    }

    /// Computes sum of 2 scalar multiplications.
    /// Faster than doing the scalar multiplications individually and then adding them. Uses lookup table
    /// returns self*a + h*b
    pub fn binary_scalar_mul(&self, h: &Self, a: &FieldElement, b: &FieldElement) -> Self {
        self.value.mul2(&a.to_bignum(), &h.to_ecp(), &b.to_bignum()).into()
    }

    /// Multiply point on the curve (element of group G1) with a scalar. Variable time operation
    /// Uses wNAF.
    pub fn scalar_mul_variable_time(&self, a: &FieldElement) -> Self {
        // TODO: Optimization: Attach the lookup table to the struct
        let table = NafLookupTable5::from(self);
        let wnaf = a.to_wnaf(5);
        GroupElement::wnaf_mul(&table, &wnaf)
    }

    pub fn double(&self) -> Self {
        let mut d = self.value.clone();
        d.dbl();
        d.into()
    }

    pub fn double_mut(&mut self) {
        self.value.dbl();
    }

    // Return multiples of itself. eg. Given `n`=5, returns self, 2*self, 3*self, 4*self, 5*self
    pub fn get_multiples(&self, n: usize) -> Vec<GroupElement> {
        let mut res = vec![self.clone()];
        for i in 2..=n {
            res.push(res[i-2] + self);
        }
        res
    }

    pub fn to_wnaf_lookup_table(&self, width: usize) -> NafLookupTable5 {
        // Only supporting table of width 5 for now
        debug_assert_eq!(width, 5);
        NafLookupTable5::from(self)
    }

    pub fn wnaf_mul(table: &NafLookupTable5, wnaf: &[i8]) -> Self {
        let mut result = GroupElement::identity();

        for n in wnaf.iter().rev() {
            result = result.double();

            let v = *n;
            if v > 0 {
                result = result + table.select(v as usize);
            } else if v < 0 {
                result = result - table.select(-v as usize);
            }
        }

        result
    }

    pub fn to_hex(&self) -> String {
        self.to_ecp().tostring()
    }
}

impl From<GroupG1> for GroupElement {
    fn from(x: GroupG1) -> Self {
        Self {
            value: x
        }
    }
}

impl From<&GroupG1> for GroupElement {
    fn from(x: &GroupG1) -> Self {
        Self {
            value: x.clone()
        }
    }
}

impl PartialEq for GroupElement {
    fn eq(&self, other: &GroupElement) -> bool {
        let mut l = self.clone();
        let mut r = other.clone();
        l.value.equals(&mut r.value)
    }
}

impl Add for GroupElement {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        self.plus(&other)
    }
}

impl Add<GroupElement> for &GroupElement {
    type Output = GroupElement;

    fn add(self, other: GroupElement) -> GroupElement {
        self.plus(&other)
    }
}

impl<'a> Add<&'a GroupElement> for GroupElement {
    type Output = Self;
    fn add(self, other: &'a GroupElement) -> Self {
        self.plus(other)
    }
}

impl AddAssign for GroupElement {
    fn add_assign(&mut self, other: Self) {
        self.add_assign_(&other)
    }
}

impl Sub for GroupElement {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        self.minus(&other)
    }
}

impl Mul<FieldElement> for GroupElement {
    type Output = Self;

    fn mul(self, other: FieldElement) -> Self {
        self.scalar_mul_const_time(&other)
    }
}

impl Mul<&FieldElement> for GroupElement {
    type Output = Self;

    fn mul(self, other: &FieldElement) -> Self {
        self.scalar_mul_const_time(other)
    }
}

impl Mul<FieldElement> for &GroupElement {
    type Output = GroupElement;

    fn mul(self, other: FieldElement) -> GroupElement {
        self.scalar_mul_const_time(&other)
    }
}

impl Mul<&FieldElement> for &GroupElement {
    type Output = GroupElement;

    fn mul(self, other: &FieldElement) -> GroupElement {
        self.scalar_mul_const_time(other)
    }
}

#[derive(Clone, Debug)]
pub struct GroupElementVector {
    elems: Vec<GroupElement>
}

impl GroupElementVector {
    pub fn new(size: usize) -> Self {
        Self {
            elems: (0..size).map(|_| GroupElement::new()).collect()
        }
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            elems: Vec::<GroupElement>::with_capacity(capacity)
        }
    }

    pub fn as_slice(&self) -> &[GroupElement] {
        &self.elems
    }

    pub fn len(&self) -> usize {
        self.elems.len()
    }

    pub fn push(&mut self, value: GroupElement) {
        self.elems.push(value)
    }

    pub fn append(&mut self, other: &mut Self) {
        self.elems.append(&mut other.elems)
    }

    /// Compute sum of all elements of a vector
    pub fn sum(&self) -> GroupElement {
        let mut accum = GroupElement::new();
        for i in 0..self.len() {
            accum += self[i];
        }
        accum
    }

    /// Multiply each field element of the vector with another given field
    /// element `n` (scale the vector)
    pub fn scale(&mut self, n: &FieldElement) {
        for i in 0..self.len() {
            self[i] = self[i] * n;
        }
    }

    pub fn scaled_by(&self, n: &FieldElement) -> Self {
        let mut scaled = Self::with_capacity(self.len());
        for i in 0..self.len() {
            scaled.push(self[i] * n)
        }
        scaled.into()
    }

    /// Computes inner product of 2 vectors, one of field elements and other of group elements.
    /// [a1, a2, a3, ...field elements].[b1, b2, b3, ...group elements] = (a1*b1 + a2*b2 + a3*b3)
    pub fn inner_product_const_time(&self, b: &FieldElementVector) -> Result<GroupElement, ValueError> {
        self.multi_scalar_mul_const_time(b)
    }

    pub fn inner_product_var_time(&self, b: &FieldElementVector) -> Result<GroupElement, ValueError> {
        self.multi_scalar_mul_var_time(b)
    }

    /// Calculates Hadamard product of 2 group element vectors.
    /// Hadamard product of `a` and `b` = `a` o `b` = (a0 o b0, a1 o b1, ...).
    /// Here `o` denotes group operation, which in elliptic curve is point addition
    pub fn hadamard_product(&self, b: &Self) -> Result<Self, ValueError> {
        check_vector_size_for_equality!(self, b)?;
        let mut hadamard_product = Self::with_capacity(self.len());
        for i in 0..self.len() {
            hadamard_product.push(self[i].plus(&b[i]));
        }
        Ok(hadamard_product)
    }

    pub fn split_at(&self, mid: usize) -> (Self, Self) {
        let (l, r) = self.as_slice().split_at(mid);
        (Self::from(l), Self::from(r))
    }

    /// Constant time multi-scalar multiplication. Naive approach computing `n` scalar multiplications and 1 addition for `n` field elements
    pub fn multi_scalar_mul_const_time_naive(&self, field_elems: &FieldElementVector) -> Result<GroupElement, ValueError> {
        check_vector_size_for_equality!(field_elems, self)?;
        let mut accum = GroupElement::new();
        for i in 0..self.len() {
            accum += self[i] * field_elems[i];
        }
        Ok(accum)
    }

    /// Constant time multi-scalar multiplication
    pub fn multi_scalar_mul_const_time(&self, field_elems: &FieldElementVector) -> Result<GroupElement, ValueError> {
        Self::_multi_scalar_mul_const_time(&self, field_elems)
    }

    /// Variable time multi-scalar multiplication
    pub fn multi_scalar_mul_var_time(&self, field_elems: &FieldElementVector) -> Result<GroupElement, ValueError> {
        Self::_multi_scalar_mul_var_time(&self, field_elems)
    }

    /// Strauss multi-scalar multiplication
    fn _multi_scalar_mul_var_time(group_elems: &GroupElementVector, field_elems: &FieldElementVector) -> Result<GroupElement, ValueError> {
        check_vector_size_for_equality!(field_elems, group_elems)?;
        let lookup_tables: Vec<_> = group_elems.as_slice()
            .into_iter()
            .map(|e| NafLookupTable5::from(e))
            .collect();

        Self::multi_scalar_mul_var_time_with_precomputation_done(&lookup_tables, field_elems)
    }

    /// Strauss multi-scalar multiplication. Passing the lookup tables since in lot of cases generators will be fixed
    pub fn multi_scalar_mul_var_time_with_precomputation_done(lookup_tables: &[NafLookupTable5],
                                                              field_elems: &FieldElementVector) -> Result<GroupElement, ValueError> {
        // Redundant check when called from multi_scalar_mul_var_time
        check_vector_size_for_equality!(field_elems, lookup_tables)?;

        let mut nafs: Vec<_> = field_elems.as_slice()
            .into_iter()
            .map(|e| e.to_wnaf(5))
            .collect();

        // Pad the NAFs with 0 so that all nafs are of same length
        let new_length = pad_collection!(nafs, 0);

        let mut r = GroupElement::identity();

        for i in (0..new_length).rev() {
            let mut t = r.double();

            for (naf, lookup_table) in nafs.iter().zip(lookup_tables.iter()) {
                if naf[i] > 0 {
                    t = t + lookup_table.select(naf[i] as usize);
                } else if naf[i] < 0 {
                    t = t - lookup_table.select(-naf[i] as usize);
                }
            }

            r = t;
        }

        Ok(r)
    }

    /// Constant time multi-scalar multiplication.
    /// Taken from Guide to Elliptic Curve Cryptography book, "Algorithm 3.48 Simultaneous multiple point multiplication" without precomputing the addition
    /// Still helps with reducing doublings
    fn _multi_scalar_mul_const_time(group_elems: &GroupElementVector, field_elems: &FieldElementVector) -> Result<GroupElement, ValueError> {
        check_vector_size_for_equality!(field_elems, group_elems)?;

        // Choosing window of size 3.
        let group_elem_multiples: Vec<_> = group_elems.as_slice()
            .into_iter()
            .map(|e| e.get_multiples(7))       // 2^3 - 1
            .collect();

        Self::multi_scalar_mul_const_time_with_precomputation_done(&group_elem_multiples, field_elems)
    }

    pub fn multi_scalar_mul_const_time_with_precomputation_done(group_elem_multiples: &Vec<Vec<GroupElement>>,
                                                              field_elems: &FieldElementVector) -> Result<GroupElement, ValueError> {
        // Redundant check when called from multi_scalar_mul_const_time
        check_vector_size_for_equality!(group_elem_multiples, field_elems)?;

        // TODO: The test shows that precomputing multiples does not help much. Experiment with bigger window.

        let mut field_elems_base_repr: Vec<_> = field_elems.as_slice()
            .into_iter()
            .map(|e| e.to_power_of_2_base(3))
            .collect();

        // Pad the representations with 0 so that all are of same length
        let new_length = pad_collection!(field_elems_base_repr, 0);

        let mut r = GroupElement::new();
        for i in (0..new_length).rev() {
            // r = r * 2^3
            r.double_mut();
            r.double_mut();
            r.double_mut();
            for (b, m) in field_elems_base_repr.iter().zip(group_elem_multiples.iter()) {
                // TODO: The following can be replaced with a pre-computation.
                if b[i] != 0 {
                    r = r + m[(b[i]-1) as usize]
                }
            }
        }
        Ok(r)
    }

    pub fn iter(&self) -> Iter<GroupElement> {
        self.as_slice().iter()
    }

}

impl From<Vec<GroupElement>> for GroupElementVector {
    fn from(x: Vec<GroupElement>) -> Self {
        Self {
            elems: x
        }
    }
}

impl From<&[GroupElement]> for GroupElementVector {
    fn from(x: &[GroupElement]) -> Self {
        Self {
            elems: x.to_vec()
        }
    }
}

impl Index<usize> for GroupElementVector {
    type Output = GroupElement;

    fn index(&self, idx: usize) -> &GroupElement {
        &self.elems[idx]
    }
}

impl IndexMut<usize> for GroupElementVector {

    fn index_mut(&mut self, idx: usize) -> &mut GroupElement {
        &mut self.elems[idx]
    }
}

impl PartialEq for GroupElementVector {
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

impl IntoIterator for GroupElementVector {
    type Item = GroupElement;
    type IntoIter = ::std::vec::IntoIter<GroupElement>;

    fn into_iter(self) -> Self::IntoIter {
        self.elems.into_iter()
    }
}

pub struct NafLookupTable5([GroupElement; 8]);

impl NafLookupTable5 {
    /// Given public A and odd x with 0 < x < 2^4, return x.A.
    pub fn select(&self, x: usize) -> GroupElement {
        debug_assert_eq!(x & 1, 1);
        debug_assert!(x < 16);

        self.0[x / 2]
    }
}

impl<'a> From<&'a GroupElement> for NafLookupTable5 {
    fn from(A: &'a GroupElement) -> Self {
        let mut Ai = [GroupElement::new(); 8];
        let A2 = A.double();
        Ai[0] = A.clone();
        for i in 0..7 {
            Ai[i + 1] = Ai[i] + A2;
        }
        // Now Ai = [A, 3A, 5A, 7A, 9A, 11A, 13A, 15A]
        Self(Ai)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::time::{Duration, Instant};

    #[test]
    fn test_scalar_mult_operators() {
        for _ in 0..10 {
            let g = GroupElement::random(None);
            let f = FieldElement::random(None);
            let m = g.scalar_mul_const_time(&f);
            // Operands can be in any order
            assert_eq!(m, g * f);
            assert_eq!(m, f * g);
        }
    }

    #[test]
    fn test_binary_scalar_mul() {
        for _ in 0..10 {
            let a = FieldElement::random(None);
            let b = FieldElement::random(None);
            let g = GroupElement::random(None);
            let h = GroupElement::random(None);
            assert_eq!(g * a + h * b, g.binary_scalar_mul(&h, &a, &b))
        }
    }

    #[test]
    fn test_group_elem_addition() {
        let a = GroupElement::random(None);
        let b = GroupElement::random(None);
        let c = GroupElement::random(None);

        let mut sum =  a + b + c;

        let mut expected_sum = GroupElement::new();
        expected_sum = expected_sum.plus(&a);
        expected_sum = expected_sum.plus(&b);
        expected_sum = expected_sum.plus(&c);
        assert_eq!(sum, expected_sum);
    }

    #[test]
    fn test_multiples() {
        for _ in 0..10 {
            let a = GroupElement::random(None);
            let mults = a.get_multiples(17);
            for i in 1..=17 {
                assert_eq!(mults[i-1], a * FieldElement::from(i as u8));
            }
        }
    }

    #[test]
    fn test_NafLookupTable5() {
        let a = GroupElement::random(None);
        let x = [1, 3, 5, 7, 9, 11, 13, 15];
        let table = NafLookupTable5::from(&a);
        for i in x.iter() {
            let f = FieldElement::from(*i as u8);
            let expected = a * f;
            assert_eq!(expected, table.select(*i as usize));
        }
    }

    #[test]
    fn test_wnaf_mul() {
        for _ in 0..100 {
            let mut a = GroupElement::random(None);
            let r = FieldElement::random(None);
            let expected = a * r;

            let table = NafLookupTable5::from(&a);
            let wnaf = r.to_wnaf(5);
            let p = GroupElement::wnaf_mul(&table, &wnaf);

            assert_eq!(expected, p);
        }
    }

    #[test]
    fn test_multi_scalar_multiplication() {
        for _ in 0..5 {
            let mut fs = vec![];
            let mut gs = vec![];
            let gen: GroupElement = GroupElement::generator();

            for i in 0..70 {
                fs.push(FieldElement::random(None));
                gs.push(gen.scalar_mul_const_time(&fs[i]));
            }

            let gv = GroupElementVector::from(gs.as_slice());
            let fv = FieldElementVector::from(fs.as_slice());
            let res = gv.multi_scalar_mul_const_time_naive(&fv).unwrap();

            let res_1 = gv.multi_scalar_mul_var_time(&fv).unwrap();

            let mut expected = GroupElement::new();
            let mut expected_1 = GroupElement::new();
            for i in 0..fs.len() {
                expected.add_assign_(&gs[i].scalar_mul_const_time(&fs[i]));
                expected_1.add_assign_(&(gs[i] * &fs[i]));
            }

            let res_2 = GroupElementVector::_multi_scalar_mul_const_time(&gv, &fv).unwrap();

            assert_eq!(expected, res);
            assert_eq!(expected_1, res);
            assert_eq!(res_1, res);
            assert_eq!(res_2, res);
        }
    }

    #[test]
    fn timing_multi_scalar_multiplication() {
        let mut fs = vec![];
        let mut gs = vec![];

        let n = 64;

        for _ in 0..n {
            fs.push(FieldElement::random(None));
            gs.push(GroupElement::random(None));
        }

        let gv = GroupElementVector::from(gs.as_slice());
        let fv = FieldElementVector::from(fs.as_slice());

        let mut start = Instant::now();
        let res = gv.multi_scalar_mul_const_time_naive(&fv).unwrap();
        let const_time_naive = start.elapsed();

        start = Instant::now();
        let res_1 = gv.multi_scalar_mul_var_time(&fv).unwrap();
        let var_time = start.elapsed();

        assert_eq!(res_1, res);

        let lookup_tables: Vec<_> = gv.as_slice()
            .into_iter()
            .map(|e| e.to_wnaf_lookup_table(5))
            .collect();

        start = Instant::now();
        let res_2 = GroupElementVector::multi_scalar_mul_var_time_with_precomputation_done(&lookup_tables, &fv).unwrap();
        let var_precomp_time = start.elapsed();

        assert_eq!(res_2, res);

        start = Instant::now();
        let res_3 = gv.multi_scalar_mul_const_time(&fv).unwrap();
        let const_time = start.elapsed();

        assert_eq!(res_3, res);

        let group_elem_multiples: Vec<_> = gv.as_slice()
            .into_iter()
            .map(|e| e.get_multiples(7))
            .collect();

        start = Instant::now();
        let res_4 = GroupElementVector::multi_scalar_mul_const_time_with_precomputation_done(&group_elem_multiples, &fv).unwrap();
        let const_precomp_time = start.elapsed();

        assert_eq!(res_4, res);

        println!("Constant time for {} size multi-scalar multiplications using naive method: {:?}", n, const_time_naive);
        println!("Constant time for {} size multi-scalar multiplications: {:?}", n, const_time);
        println!("Constant time with pre-computation for {} size multi-scalar multiplications: {:?}", n, const_precomp_time);
        println!("Variable time for {} size multi-scalar multiplications: {:?}", n, var_time);
        println!("Variable time with pre-computation for {} size multi-scalar multiplications: {:?}", n, var_precomp_time);
    }

    #[test]
    fn timing_wnaf_mul() {
        let mut fs = vec![];
        let mut gs = vec![];

        let n = 64;
        let w = 5;

        for _ in 0..n {
            fs.push(FieldElement::random(None));
            gs.push(GroupElement::random(None));
        }

        let gv = GroupElementVector::from(gs.as_slice());
        let fv = FieldElementVector::from(fs.as_slice());

        let mut start = Instant::now();
        for i in 0..n {
            // The compiler might not execute the statement below
            gv[i] * fv[i];
        }
        println!("Time for {} scalar multiplications: {:?}", n, start.elapsed());

        start = Instant::now();
        for i in 0..n {
            let naf = fv[i].to_wnaf(w);
            let table = gv[i].to_wnaf_lookup_table(w);
            // The compiler might not execute the statement below
            GroupElement::wnaf_mul(&table, &naf);
        }
        println!("Time for {} scalar multiplications using wnaf: {:?}", n, start.elapsed());
    }

    #[test]
    fn timing_group_elem_addition() {
        let count = 100;
        let points: Vec<GroupElement> = (0..100).map(|_| GroupElement::random(None)).collect();
        let mut R = GroupElement::random(None);
        let mut start = Instant::now();
        for i in 0..count {
            R = R + points[i];
        }
        println!("Addition time for {} elems = {:?}", count, start.elapsed());
    }
}
