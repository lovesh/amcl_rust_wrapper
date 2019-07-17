use crate::constants::{CurveOrder, GroupG1_SIZE};
use crate::errors::{SerzDeserzError, ValueError};
use crate::field_elem::{FieldElement, FieldElementVector};
use crate::group_elem::{GroupElement, GroupElementVector};
use crate::types::GroupG1;
use crate::utils::hash_msg;
use std::ops::{Add, AddAssign, Index, IndexMut, Mul, Neg, Sub};

use std::fmt;
use std::hash::{Hash, Hasher};
use std::slice::Iter;

#[derive(Clone, Debug)]
pub struct G1 {
    value: GroupG1,
}

impl GroupElement for G1 {
    fn new() -> Self {
        Self {
            value: GroupG1::new(),
        }
    }

    fn identity() -> Self {
        let mut v = GroupG1::new();
        v.inf();
        Self { value: v }
    }

    fn generator() -> Self {
        GroupG1::generator().into()
    }

    fn is_identity(&self) -> bool {
        self.value.is_infinity()
    }

    fn set_to_identity(&mut self) {
        self.value.inf()
    }

    fn from_msg_hash(msg: &[u8]) -> Self {
        GroupG1::mapit(&hash_msg(msg)).into()
    }

    fn to_bytes(&self) -> Vec<u8> {
        let mut temp = GroupG1::new();
        temp.copy(&self.value);
        let mut bytes: [u8; GroupG1_SIZE] = [0; GroupG1_SIZE];
        temp.tobytes(&mut bytes, false);
        bytes.to_vec()
    }

    fn from_bytes(bytes: &[u8]) -> Result<Self, SerzDeserzError> {
        if bytes.len() != GroupG1_SIZE {
            return Err(SerzDeserzError::G1BytesIncorrectSize(
                bytes.len(),
                GroupG1_SIZE,
            ));
        }
        Ok(GroupG1::frombytes(bytes).into())
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
        return false;
    }

    fn has_correct_order(&self) -> bool {
        return self.value.mul(&CurveOrder).is_infinity();
    }
}

/// Represents an element of the sub-group of the elliptic curve over the prime field
impl G1 {
    /// Return underlying elliptic curve point, ECP
    pub fn to_ecp(&self) -> GroupG1 {
        self.value.clone()
    }

    /// Computes sum of 2 scalar multiplications.
    /// Faster than doing the scalar multiplications individually and then adding them. Uses lookup table
    /// returns self*a + h*b
    pub fn binary_scalar_mul(&self, h: &Self, a: &FieldElement, b: &FieldElement) -> Self {
        self.value
            .mul2(&a.to_bignum(), &h.to_ecp(), &b.to_bignum())
            .into()
    }

    /// Multiply point on the curve (element of group G1) with a scalar. Variable time operation
    /// Uses wNAF.
    pub fn scalar_mul_variable_time(&self, a: &FieldElement) -> Self {
        // TODO: Optimization: Attach the lookup table to the struct
        let table = G1LookupTable::from(self);
        let wnaf = a.to_wnaf(5);
        G1::wnaf_mul(&table, &wnaf)
    }

    /// Return multiples of itself. eg. Given `n`=5, returns self, 2*self, 3*self, 4*self, 5*self
    pub fn get_multiples(&self, n: usize) -> Vec<G1> {
        // TODO: Can use `selector` from ECP
        let mut res = vec![self.clone()];
        for i in 2..=n {
            res.push(&res[i - 2] + self);
        }
        res
    }

    pub fn to_wnaf_lookup_table(&self, width: usize) -> G1LookupTable {
        // Only supporting table of width 5 for now
        debug_assert_eq!(width, 5);
        G1LookupTable::from(self)
    }

    pub fn wnaf_mul(table: &G1LookupTable, wnaf: &[i8]) -> Self {
        let mut result = G1::identity();

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
}

impl_group_elem_traits!(G1);

impl_group_elem_conversions!(G1, GroupG1, GroupG1_SIZE);

impl_group_elem_ops!(G1);

impl_scalar_mul_ops!(G1);

#[derive(Clone, Debug)]
pub struct G1Vector {
    elems: Vec<G1>,
}

impl_group_elem_vec_ops!(G1, G1Vector);

impl_group_elem_vec_conversions!(G1, G1Vector);

impl G1Vector {
    /// Computes inner product of 2 vectors, one of field elements and other of group elements.
    /// [a1, a2, a3, ...field elements].[b1, b2, b3, ...group elements] = (a1*b1 + a2*b2 + a3*b3)
    pub fn inner_product_const_time(&self, b: &FieldElementVector) -> Result<G1, ValueError> {
        self.multi_scalar_mul_const_time(b)
    }

    pub fn inner_product_var_time(&self, b: &FieldElementVector) -> Result<G1, ValueError> {
        self.multi_scalar_mul_var_time(b)
    }

    /// Calculates Hadamard product of 2 group element vectors.
    /// Hadamard product of `a` and `b` = `a` o `b` = (a0 o b0, a1 o b1, ...).
    /// Here `o` denotes group operation, which in elliptic curve is point addition
    pub fn hadamard_product(&self, b: &Self) -> Result<Self, ValueError> {
        check_vector_size_for_equality!(self, b)?;
        let mut hadamard_product = Self::with_capacity(self.len());
        for i in 0..self.len() {
            hadamard_product.push(&self[i] + &b[i]);
        }
        Ok(hadamard_product)
    }

    pub fn split_at(&self, mid: usize) -> (Self, Self) {
        let (l, r) = self.as_slice().split_at(mid);
        (Self::from(l), Self::from(r))
    }

    /// Constant time multi-scalar multiplication. Naive approach computing `n` scalar multiplications and 1 addition for `n` field elements
    pub fn multi_scalar_mul_const_time_naive(
        &self,
        field_elems: &FieldElementVector,
    ) -> Result<G1, ValueError> {
        check_vector_size_for_equality!(field_elems, self)?;
        let mut accum = G1::new();
        for i in 0..self.len() {
            accum += &self[i] * &field_elems[i];
        }
        Ok(accum)
    }

    /// Constant time multi-scalar multiplication
    pub fn multi_scalar_mul_const_time(
        &self,
        field_elems: &FieldElementVector,
    ) -> Result<G1, ValueError> {
        Self::_multi_scalar_mul_const_time(&self, field_elems)
    }

    /// Variable time multi-scalar multiplication
    pub fn multi_scalar_mul_var_time(
        &self,
        field_elems: &FieldElementVector,
    ) -> Result<G1, ValueError> {
        Self::_multi_scalar_mul_var_time(&self, field_elems)
    }

    /// Strauss multi-scalar multiplication
    fn _multi_scalar_mul_var_time(
        group_elems: &G1Vector,
        field_elems: &FieldElementVector,
    ) -> Result<G1, ValueError> {
        check_vector_size_for_equality!(field_elems, group_elems)?;
        let lookup_tables: Vec<_> = group_elems
            .as_slice()
            .into_iter()
            .map(|e| G1LookupTable::from(e))
            .collect();

        Self::multi_scalar_mul_var_time_with_precomputation_done(&lookup_tables, field_elems)
    }

    /// Strauss multi-scalar multiplication. Passing the lookup tables since in lot of cases generators will be fixed
    pub fn multi_scalar_mul_var_time_with_precomputation_done(
        lookup_tables: &[G1LookupTable],
        field_elems: &FieldElementVector,
    ) -> Result<G1, ValueError> {
        // Redundant check when called from multi_scalar_mul_var_time
        check_vector_size_for_equality!(field_elems, lookup_tables)?;

        let mut nafs: Vec<_> = field_elems
            .as_slice()
            .into_iter()
            .map(|e| e.to_wnaf(5))
            .collect();

        // Pad the NAFs with 0 so that all nafs are of same length
        let new_length = pad_collection!(nafs, 0);

        let mut r = G1::identity();

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
    fn _multi_scalar_mul_const_time(
        group_elems: &G1Vector,
        field_elems: &FieldElementVector,
    ) -> Result<G1, ValueError> {
        check_vector_size_for_equality!(field_elems, group_elems)?;

        // Choosing window of size 3.
        let group_elem_multiples: Vec<_> = group_elems
            .as_slice()
            .into_iter()
            .map(|e| e.get_multiples(7)) // 2^3 - 1
            .collect();

        Self::multi_scalar_mul_const_time_with_precomputation_done(
            &group_elem_multiples,
            field_elems,
        )
    }

    pub fn multi_scalar_mul_const_time_with_precomputation_done(
        group_elem_multiples: &Vec<Vec<G1>>,
        field_elems: &FieldElementVector,
    ) -> Result<G1, ValueError> {
        // Redundant check when called from multi_scalar_mul_const_time
        check_vector_size_for_equality!(group_elem_multiples, field_elems)?;

        // TODO: The test shows that precomputing multiples does not help much. Experiment with bigger window.

        let mut field_elems_base_repr: Vec<_> = field_elems
            .as_slice()
            .into_iter()
            .map(|e| e.to_power_of_2_base(3))
            .collect();

        // Pad the representations with 0 so that all are of same length
        let new_length = pad_collection!(field_elems_base_repr, 0);

        let mut r = G1::new();
        for i in (0..new_length).rev() {
            // r = r * 2^3
            r.double_mut();
            r.double_mut();
            r.double_mut();
            for (b, m) in field_elems_base_repr
                .iter()
                .zip(group_elem_multiples.iter())
            {
                // TODO: The following can be replaced with a pre-computation.
                if b[i] != 0 {
                    r = r + &m[(b[i] - 1) as usize]
                }
            }
        }
        Ok(r)
    }
}

pub struct G1LookupTable([G1; 8]);

impl G1LookupTable {
    /// Given public A and odd x with 0 < x < 2^4, return x.A.
    pub fn select(&self, x: usize) -> &G1 {
        debug_assert_eq!(x & 1, 1);
        debug_assert!(x < 16);

        &self.0[x / 2]
    }
}

impl<'a> From<&'a G1> for G1LookupTable {
    fn from(A: &'a G1) -> Self {
        let mut Ai: [G1; 8] = [G1::new(), G1::new(), G1::new(), G1::new(), G1::new(), G1::new(), G1::new(), G1::new()];
        let A2 = A.double();
        Ai[0] = A.clone();
        for i in 0..7 {
            Ai[i + 1] = &Ai[i] + &A2;
        }
        // Now Ai = [A, 3A, 5A, 7A, 9A, 11A, 13A, 15A]
        Self(Ai)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::collections::{HashMap, HashSet};
    use std::time::{Duration, Instant};

    #[test]
    fn test_binary_scalar_mul() {
        for _ in 0..10 {
            let a = FieldElement::random();
            let b = FieldElement::random();
            let g = G1::random();
            let h = G1::random();
            assert_eq!(&g * &a + &h * &b, g.binary_scalar_mul(&h, &a, &b))
        }
    }

    #[test]
    fn test_multiples() {
        for _ in 0..10 {
            let a = G1::random();
            let mults = a.get_multiples(17);
            for i in 1..=17 {
                assert_eq!(mults[i - 1], (&a * FieldElement::from(i as u8)));
            }
        }
    }

    #[test]
    fn test_NafLookupTable5() {
        let a = G1::random();
        let x = [1, 3, 5, 7, 9, 11, 13, 15];
        let table = G1LookupTable::from(&a);
        for i in x.iter() {
            let f = FieldElement::from(*i as u8);
            let expected = &a * f;
            assert_eq!(expected, *table.select(*i as usize));
        }
    }

    #[test]
    fn test_wnaf_mul() {
        for _ in 0..100 {
            let a = G1::random();
            let r = FieldElement::random();
            let expected = &a * &r;

            let table = G1LookupTable::from(&a);
            let wnaf = r.to_wnaf(5);
            let p = G1::wnaf_mul(&table, &wnaf);

            assert_eq!(expected, p);
        }
    }

    #[test]
    fn test_multi_scalar_multiplication() {
        for _ in 0..5 {
            let mut fs = vec![];
            let mut gs = vec![];
            let gen: G1 = G1::generator();

            for i in 0..70 {
                fs.push(FieldElement::random());
                gs.push(gen.scalar_mul_const_time(&fs[i]));
            }

            let gv = G1Vector::from(gs.as_slice());
            let fv = FieldElementVector::from(fs.as_slice());
            let res = gv.multi_scalar_mul_const_time_naive(&fv).unwrap();

            let res_1 = gv.multi_scalar_mul_var_time(&fv).unwrap();

            let mut expected = G1::new();
            let mut expected_1 = G1::new();
            for i in 0..fs.len() {
                expected.add_assign_(&gs[i].scalar_mul_const_time(&fs[i]));
                expected_1.add_assign_(&(&gs[i] * &fs[i]));
            }

            let res_2 = G1Vector::_multi_scalar_mul_const_time(&gv, &fv).unwrap();

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
            fs.push(FieldElement::random());
            gs.push(G1::random());
        }

        let gv = G1Vector::from(gs.as_slice());
        let fv = FieldElementVector::from(fs.as_slice());

        let mut start = Instant::now();
        let res = gv.multi_scalar_mul_const_time_naive(&fv).unwrap();
        let const_time_naive = start.elapsed();

        start = Instant::now();
        let res_1 = gv.multi_scalar_mul_var_time(&fv).unwrap();
        let var_time = start.elapsed();

        assert_eq!(res_1, res);

        let lookup_tables: Vec<_> = gv
            .as_slice()
            .into_iter()
            .map(|e| e.to_wnaf_lookup_table(5))
            .collect();

        start = Instant::now();
        let res_2 =
            G1Vector::multi_scalar_mul_var_time_with_precomputation_done(&lookup_tables, &fv)
                .unwrap();
        let var_precomp_time = start.elapsed();

        assert_eq!(res_2, res);

        start = Instant::now();
        let res_3 = gv.multi_scalar_mul_const_time(&fv).unwrap();
        let const_time = start.elapsed();

        assert_eq!(res_3, res);

        let group_elem_multiples: Vec<_> = gv
            .as_slice()
            .into_iter()
            .map(|e| e.get_multiples(7))
            .collect();

        start = Instant::now();
        let res_4 = G1Vector::multi_scalar_mul_const_time_with_precomputation_done(
            &group_elem_multiples,
            &fv,
        )
        .unwrap();
        let const_precomp_time = start.elapsed();

        assert_eq!(res_4, res);

        println!(
            "Constant time for {} size multi-scalar multiplications using naive method: {:?}",
            n, const_time_naive
        );
        println!(
            "Constant time for {} size multi-scalar multiplications: {:?}",
            n, const_time
        );
        println!(
            "Constant time with pre-computation for {} size multi-scalar multiplications: {:?}",
            n, const_precomp_time
        );
        println!(
            "Variable time for {} size multi-scalar multiplications: {:?}",
            n, var_time
        );
        println!(
            "Variable time with pre-computation for {} size multi-scalar multiplications: {:?}",
            n, var_precomp_time
        );
    }

    #[test]
    fn timing_wnaf_mul() {
        let mut fs = vec![];
        let mut gs = vec![];

        let n = 64;
        let w = 5;

        for _ in 0..n {
            fs.push(FieldElement::random());
            gs.push(G1::random());
        }

        let gv = G1Vector::from(gs.as_slice());
        let fv = FieldElementVector::from(fs.as_slice());

        let mut start = Instant::now();
        for i in 0..n {
            // The compiler might not execute the statement below
            let _ = &gv[i] * &fv[i];
        }
        println!(
            "Time for {} scalar multiplications: {:?}",
            n,
            start.elapsed()
        );

        start = Instant::now();
        for i in 0..n {
            let naf = fv[i].to_wnaf(w);
            let table = gv[i].to_wnaf_lookup_table(w);
            // The compiler might not execute the statement below
            G1::wnaf_mul(&table, &naf);
        }
        println!(
            "Time for {} scalar multiplications using wnaf: {:?}",
            n,
            start.elapsed()
        );
    }
}
