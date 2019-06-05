#![allow(non_snake_case)]

extern crate amcl;

#[macro_use]
extern crate lazy_static;

pub use amcl::bls381 as ECCurve;
//pub use amcl::bn254 as ECCurve;


pub mod types;
pub mod constants;

#[macro_use]
pub mod errors;

#[macro_use]
pub mod macros;

pub mod utils;

#[macro_use]
pub mod field_elem;
#[macro_use]
pub mod group_elem;
pub mod group_elem_g1;
pub mod commitment;