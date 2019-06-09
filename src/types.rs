use super::ECCurve::big::BIG;
use super::ECCurve::dbig::DBIG;
use super::ECCurve::ecp::ECP;
use super::ECCurve::ecp2::ECP2;
use super::ECCurve::fp12::FP12;

pub type BigNum = BIG;
pub type DoubleBigNum = DBIG;
pub type GroupG1 = ECP;
pub type GroupG2 = ECP2;
// TODO: Implement GroupElement for Gt
pub type GroupGT = FP12;
