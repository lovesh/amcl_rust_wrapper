#[allow(non_snake_case)]
#[allow(non_upper_case_globals)]

pub extern crate amcl;

#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate arrayref;

#[cfg(feature = "bn254")]
pub use amcl::bn254 as ECCurve;

#[cfg(feature = "bls381")]
pub use amcl::bls381 as ECCurve;

#[cfg(feature = "secp256k1")]
pub use amcl::secp256k1 as ECCurve;

#[cfg(feature = "ed25519")]
pub use amcl::ed25519 as ECCurve;

pub mod constants;
pub mod types;

#[macro_use]
pub mod errors;

#[macro_use]
pub mod macros;

pub mod utils;

#[macro_use]
pub mod curve_order_elem;
#[macro_use]
pub mod group_elem;
#[macro_use]
pub mod group_elem_g1;
pub mod commitment;
#[macro_use]
pub mod univar_poly;

#[cfg(any(feature = "bls381", feature = "bn254"))]
pub mod types_g2;

#[cfg(any(feature = "bls381", feature = "bn254"))]
#[macro_use]
pub mod group_elem_g2;

#[cfg(any(feature = "bls381", feature = "bn254"))]
#[macro_use]
pub mod extension_field_gt;

// TODO: Move the timing tests to benchmark
