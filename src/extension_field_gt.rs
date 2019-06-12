use crate::types::GroupGT;

use super::ECCurve::pair::{ate, fexp};
use super::ECCurve::fp12::FP12;
use crate::field_elem::FieldElement;
use crate::group_elem::GroupElement;
use crate::group_elem_g1::G1;
use crate::group_elem_g2::G2;
use std::fmt;


pub struct GT {
    value: GroupGT
}

impl fmt::Display for GT {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut c = self.value.clone();
        write!(f, "{}", c.tostring())
    }
}

impl GT {
    pub fn ate_pairing(g1: &G1, g2: &G2) -> Self {
        let e = ate(&g2.to_ecp(), &g1.to_ecp());
        Self { value: fexp(&e) }
    }

    pub fn mul(a: &Self, b: &Self) -> Self {
        let mut m = FP12::new_copy(&a.value);
        m.mul(&b.value);
        Self { value: m }
    }

    pub fn pow(&self, e: &FieldElement) -> Self {
        Self { value: self.value.pow(&e.to_bignum()) }
    }
}

impl PartialEq for GT {
    fn eq(&self, other: &GT) -> bool {
        self.value.equals(&other.value)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_ate_pairing() {
        let g1 = G1::random();
        let h1 = G1::random();
        let g2 = G2::random();
        let h2 = G2::random();

        // e(g1 + h1, g2) == e(g1, g2)*e(h1, g2)
        let lhs = GT::ate_pairing(&(g1 + h1), &g2);
        let rhs = GT::mul(&GT::ate_pairing(&g1, &g2), &GT::ate_pairing(&h1, &g2));
        assert!(lhs == rhs);

        // e(g1, g2+h2) == e(g1, g2)*e(g1, h2)
        let lhs = GT::ate_pairing(&g1, &(g2 + h2));
        let rhs = GT::mul(&GT::ate_pairing(&g1, &g2), &GT::ate_pairing(&g1, &h2));
        assert!(lhs == rhs);

        let r = FieldElement::random();
        // e(g1, g2^r) == e(g1^r, g2) == e(g1, g2)^r
        let p1 = GT::ate_pairing(&g1, &(g2 * r));
        let p2 = GT::ate_pairing(&(g1 * r), &g2);
        let mut p = GT::ate_pairing(&g1, &g2);
        p = p.pow(&r);
        assert!(p1 == p2);
        assert!(p1 == p);
    }
}