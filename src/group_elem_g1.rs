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

impl_group_elem_traits!(G1);

impl_group_elem_conversions!(G1, GroupG1, GroupG1_SIZE);

impl_group_elem_ops!(G1);

impl_scalar_mul_ops!(G1);

impl_group_element_lookup_table!(G1, G1LookupTable);

/// Represents an element of the sub-group of the elliptic curve over the prime field
impl_optmz_scalar_mul_ops!(G1, GroupG1, G1LookupTable);

#[derive(Clone, Debug)]
pub struct G1Vector {
    elems: Vec<G1>,
}

impl_group_elem_vec_ops!(G1, G1Vector);

impl_group_elem_vec_product_ops!(G1, G1Vector, G1LookupTable);

impl_group_elem_vec_conversions!(G1, G1Vector);

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

        let f_refs: Vec<&FieldElement> = fs.iter().map(|f| f).collect();
        start = Instant::now();
        let res_2 =
            G1Vector::multi_scalar_mul_var_time_with_precomputation_done(&lookup_tables, f_refs)
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
