#![no_main]
#[macro_use] extern crate libfuzzer_sys;
extern crate amcl_wrapper;

use amcl_wrapper::field_elem::FieldElement;
use amcl_wrapper::group_elem::GroupElement;
use amcl_wrapper::group_elem_g1::G1;
use amcl_wrapper::group_elem_g2::G2;

fuzz_target!(|data: &[u8]| {
    // fuzzed code goes here
    let _ = G1::from_bytes(data);
    let _ = G2::from_bytes(data);
    let _ = G1::from_msg_hash(data);
    let _ = G2::from_msg_hash(data);
    let _ = FieldElement::from_msg_hash(data);
});
