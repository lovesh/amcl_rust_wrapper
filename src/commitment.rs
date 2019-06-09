use crate::errors::ValueError;
use crate::field_elem::{FieldElement, FieldElementVector};
use crate::group_elem::GroupElementVector;
use crate::group_elem_g1::{G1Vector, G1};

/// Commit to field element `elem` with randomness `r` given groups elements `g` and `h`, i.e. compute g^elem.h^r
pub fn commit_to_field_element(g: &G1, h: &G1, elem: &FieldElement, r: &FieldElement) -> G1 {
    g.binary_scalar_mul(h, elem, r)
}

/// Commit to field element vectors `a` and `b` with random field element `c`
/// Given group element vectors `g` and `h` and group element `u`, compute
/// (a1*g1 + a2*g2 + a3*g3) + (b1*h1 + b2*h2 + b3*h3) + c*u
pub fn commit_to_field_element_vectors(
    g: &G1Vector,
    h: &G1Vector,
    u: &G1,
    a: &FieldElementVector,
    b: &FieldElementVector,
    c: &FieldElement,
) -> Result<G1, ValueError> {
    /*let a_g = g.inner_product_const_time(a)?;
    let b_h = h.inner_product_const_time(b)?;
    let c_u = u * c;
    Ok(a_g + b_h + c_u)*/
    let mut combined_g: Vec<G1> = vec![];
    combined_g.extend_from_slice(g.as_slice());
    combined_g.extend_from_slice(h.as_slice());
    combined_g.push(*u);

    let mut combined_f: Vec<FieldElement> = vec![];
    combined_f.extend_from_slice(a.as_slice());
    combined_f.extend_from_slice(b.as_slice());
    combined_f.push(*c);

    G1Vector::from(combined_g).inner_product_const_time(&FieldElementVector::from(combined_f))
}
