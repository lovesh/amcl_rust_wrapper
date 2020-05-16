use crate::amcl::hmac;
use crate::constants::{
    CurveOrder, GroupG1_SIZE, CURVETYPE, G1_COMP_BYTE_SIZE, HASH_TYPE, MONTGOMERY,
};
use crate::errors::{SerzDeserzError, ValueError};
use crate::field_elem::{FieldElement, FieldElementVector};
use crate::group_elem::{GroupElement, GroupElementVector};
use crate::types::{BigNum, GroupG1, FP};
use crate::utils::{hash_msg, hash_to_field};
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
pub struct G1 {
    value: GroupG1,
}

impl fmt::Debug for G1 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "ECP: [ {} ]", self.value.tostring())
    }
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

    /// This is an arbitrary choice. Any group element can be a generator
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

    impl_group_elem_byte_conversion_methods!(
        GroupG1,
        GroupG1_SIZE,
        G1_COMP_BYTE_SIZE,
        SerzDeserzError::G1BytesIncorrectSize
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
    /// Returns `(x)` or `(x,y)` depending on the curve being a Montgomery curve or not.
    /// `x` and `y` are hex representations of FP
    fn to_hex(&self) -> String {
        self.value.tostring()
    }

    fn from_hex(mut string: String) -> Result<Self, SerzDeserzError> {
        if &string == "infinity" {
            return Ok(Self::new());
        }

        unbound_bounded_string!(string, '(', ')', SerzDeserzError::CannotParseG1);

        if CURVETYPE == MONTGOMERY {
            // Only x coordinate is needed for Montgomery curves
            let x_big = FieldElement::parse_hex_as_bignum(string)?;
            Ok(Self {
                value: GroupG1::new_big(&x_big),
            })
        } else {
            // Assuming string as "x,y"
            let (x, y) = split_string_to_2_tuple!(string, SerzDeserzError::CannotParseG1);

            let x_big = FieldElement::parse_hex_as_bignum(x)?;
            let y_big = FieldElement::parse_hex_as_bignum(y)?;
            Ok(Self {
                value: GroupG1::new_bigs(&x_big, &y_big),
            })
        }
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

impl G1 {
    /// Computes sum of 2 scalar multiplications.
    /// Faster than doing the scalar multiplications individually and then adding them. Uses lookup table
    /// returns self*a + h*b
    pub fn binary_scalar_mul(&self, h: &Self, a: &FieldElement, b: &FieldElement) -> Self {
        self.value
            .mul2(&a.to_bignum(), &h.to_ecp(), &b.to_bignum())
            .into()
    }

    /// Hashes a byte slice to a group element according to the hash to curve point IETF standard
    /// https://datatracker.ietf.org/doc/draft-irtf-cfrg-hash-to-curve/?include_text=1
    /// `domain_separation_tag` should be unique between protocols as well as curves, eg. protocol A and
    /// protocol B should use different `domain_separation_tag` while hashing to the same curve and
    /// protocol A should use different `domain_separation_tag` while hashing to different curves.
    /// Look at section 3.1 of the standard for more details
    pub fn hash_to_curve(domain_separation_tag: &[u8], msg: &[u8]) -> G1 {
        // Get 2 field elements as FP
        let mut u: [FP; 2] = [FP::new(), FP::new()];
        hash_to_field(
            hmac::MC_SHA2,
            HASH_TYPE,
            domain_separation_tag,
            msg,
            &mut u,
            2,
        );

        // Map each FP to a curve point and add the points
        let mut P = GroupG1::map2point(&u[0]);
        let P1 = GroupG1::map2point(&u[1]);
        P.add(&P1);
        // clear the cofactor of the addition point
        P.cfp();

        Self { value: P }
    }
}

impl_group_elem_traits!(G1, GroupG1);

impl_group_elem_serz!(G1, GroupG1, "G1");

impl_group_elem_conversions!(G1, GroupG1, GroupG1_SIZE, G1_COMP_BYTE_SIZE);

impl_group_elem_ops!(G1);

impl_scalar_mul_ops!(G1);

impl_group_element_lookup_table!(G1, G1LookupTable);

/// Represents an element of the sub-group of the elliptic curve over the prime field
impl_optmz_scalar_mul_ops!(G1, GroupG1, G1LookupTable);

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct G1Vector {
    elems: Vec<G1>,
}

impl_group_elem_vec_ops!(G1, G1Vector);

impl_group_elem_vec_product_ops!(G1, G1Vector, G1LookupTable);

impl_group_elem_vec_conversions!(G1, G1Vector);

#[cfg(test)]
mod test {
    use super::*;
    use std::borrow::Borrow;
    use std::collections::{HashMap, HashSet};
    use std::time::{Duration, Instant};

    #[test]
    fn test_parse_hex_for_FP() {
        // TODO:
    }

    #[test]
    fn test_parse_bad_hex_for_FP() {
        // TODO:
    }

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
        let res_1 = gv.multi_scalar_mul_var_time(fv.iter()).unwrap();
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
        let res_3 = gv.multi_scalar_mul_const_time(fv.as_ref()).unwrap();
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
            fv.as_slice(),
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
