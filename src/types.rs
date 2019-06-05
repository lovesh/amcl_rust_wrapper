use super::ECCurve::big::BIG;
use super::ECCurve::dbig::DBIG;
use super::ECCurve::ecp::ECP;
use super::ECCurve::ecp2::ECP2;

pub type BigNum = BIG;
pub type DoubleBigNum = DBIG;
pub type GroupG1 = ECP;
pub type GroupG2 = ECP2;
