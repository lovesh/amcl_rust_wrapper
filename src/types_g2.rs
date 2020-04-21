use crate::constants::MODBYTES;

use super::types::GroupG2;
// Byte size of element in group G2
pub const GROUP_G2_SIZE: usize = (4 * MODBYTES) as usize;

// Byte size of element in group GT
pub const GROUP_GT_SIZE: usize = (12 * MODBYTES) as usize;

lazy_static! {
    pub static ref GENERATOR_G2: GroupG2 = GroupG2::generator();
}
