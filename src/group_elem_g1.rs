use crate::constants::{CURVE_ORDER, GROUP_G1_SIZE};
use crate::errors::{SerzDeserzError, ValueError};
use crate::field_elem::{FieldElement, FieldElementVector};
use crate::group_elem::{GroupElement, GroupElementVector};
use crate::types::{GroupG1, FP};
use crate::utils::hash_msg;
use std::ops::{Add, AddAssign, Index, IndexMut, Mul, Neg, Sub, SubAssign};

use std::fmt;
use std::hash::{Hash, Hasher};
use std::slice::Iter;

use crate::rayon::iter::IntoParallelRefMutIterator;
use rayon::prelude::*;
use serde::de::{Deserialize, Deserializer, Error as DError, Visitor};
use serde::ser::{Serialize, Serializer};
use std::str::{FromStr, SplitWhitespace};
use zeroize::Zeroize;
use hash2curve::HashToCurveXmd;

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

    fn hash_to_curve(msg: &[u8], dst: &hash2curve::DomainSeparationTag) -> Self {
        let hasher = hash2curve::bls381g1::Bls12381G1Sswu::new(dst.clone());
        match hasher.hash_to_curve_xmd::<sha3::Sha3_256, &[u8]>(msg) {
            Ok(p) => p.into(),
            Err(_) => Self::identity()
        }
    }

    fn to_bytes(&self) -> Vec<u8> {
        let mut bytes: [u8; GROUP_G1_SIZE] = [0; GROUP_G1_SIZE];
        self.write_to_slice_unchecked(&mut bytes);
        bytes.to_vec()
    }

    fn from_bytes(bytes: &[u8]) -> Result<Self, SerzDeserzError> {
        if bytes.len() != GROUP_G1_SIZE {
            return Err(SerzDeserzError::G1BytesIncorrectSize(
                bytes.len(),
                GROUP_G1_SIZE,
            ));
        }
        Ok(GroupG1::frombytes(bytes).into())
    }

    fn write_to_slice(&self, target: &mut [u8]) -> Result<(), SerzDeserzError> {
        if target.len() != GROUP_G1_SIZE {
            return Err(SerzDeserzError::G1BytesIncorrectSize(
                target.len(),
                GROUP_G1_SIZE,
            ));
        }
        self.write_to_slice_unchecked(target);
        Ok(())
    }

    fn write_to_slice_unchecked(&self, target: &mut [u8]) {
        let mut temp = GroupG1::new();
        temp.copy(&self.value);
        temp.tobytes(target, false);
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
        let x = parse_hex_as_fp(&mut iter)?;
        let y = parse_hex_as_fp(&mut iter)?;
        let z = parse_hex_as_fp(&mut iter)?;
        let mut value = GroupG1::new();
        value.setpx(x);
        value.setpy(y);
        value.setpz(z);
        Ok(Self { value })
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
        return self.value.mul(&CURVE_ORDER).is_infinity();
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
}

impl_group_elem_traits!(G1, GroupG1);

impl_group_elem_conversions!(G1, GroupG1, GROUP_G1_SIZE);

impl_group_elem_ops!(G1);

impl_scalar_mul_ops!(G1);

impl_group_element_lookup_table!(G1, G1LookupTable);

// Represents an element of the sub-group of the elliptic curve over the prime field
impl_optmz_scalar_mul_ops!(G1, GroupG1, G1LookupTable);

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct G1Vector {
    elems: Vec<G1>,
}

impl_group_elem_vec_ops!(G1, G1Vector);

impl_group_elem_vec_product_ops!(G1, G1Vector, G1LookupTable);

impl_group_elem_vec_conversions!(G1, G1Vector);

/// Parse given hex string as FP
pub fn parse_hex_as_fp(iter: &mut SplitWhitespace) -> Result<FP, SerzDeserzError> {
    // Logic almost copied from AMCL but with error handling and constant time execution.
    // Constant time is important as hex is used during serialization and deserialization.
    // A seemingly effortless solution is to filter string for errors and pad with 0s before
    // passing to AMCL but that would be expensive as the string is scanned twice
    let xes = match iter.next() {
        Some(i) => {
            // Parsing as u32 since xes cannot be negative
            match u32::from_str(i) {
                Ok(xes) => xes as i32,
                Err(_) => return Err(SerzDeserzError::CannotParseFP),
            }
        }
        None => return Err(SerzDeserzError::CannotParseFP),
    };

    let x = match iter.next() {
        Some(i) => FieldElement::parse_hex_as_bignum(i.to_string())?,
        None => return Err(SerzDeserzError::CannotParseFP),
    };

    Ok(FP { x, xes })
}

#[cfg(test)]
mod test {
    use super::*;
    use std::time::Instant;
    use hash2curve::DomainSeparationTag;

    #[test]
    fn test_hash_to_curve() {
        let e = G1::from_hex("1 0546197CDBA187E858730894C66FAEB35E7DBE4C61646786FB85B3EBB78377B1711797A884CBE8302A23463FFFD00190 1 11CF30309EA1AF4BB47FC4D5219529347F9576201EE34DE933C96F83FBB8B2AC22387B593C5F148924B571FE605B337F 2 13317C30F3A0D636D56A23C34FDD80B891ECBDE7C2B7D6E16B0F4B0B7E6D26CB6147ACDE629C4A23C57400D203A9FB84".to_string()).unwrap();
        let dst = DomainSeparationTag::new("hash_to_curve_", Some("test"), None, None).unwrap();
        let g = G1::hash_to_curve(b"message to be hashed", &dst);
        assert!(!g.is_identity());
        assert_eq!(e, g);
    }

    #[test]
    fn test_parse_hex_for_fp() {
        // TODO:
    }

    #[test]
    fn test_parse_bad_hex_for_fp() {
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
