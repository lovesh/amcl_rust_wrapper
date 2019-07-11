use crate::constants::MODBYTES;

use super::types::GroupG2;
// Byte size of element in group G2
pub const GroupG2_SIZE: usize = (4 * MODBYTES) as usize;

lazy_static! {
    pub static ref GeneratorG2: GroupG2 = GroupG2::generator();
}
