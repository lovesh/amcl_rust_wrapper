use crate::constants::MODBYTES;

use super::types::GroupG2;

// Byte size of uncompressed element in group G2, 1 extra byte for compression flag
pub const GroupG2_SIZE: usize = (4 * MODBYTES + 1) as usize;

// Byte size of compressed element in group G2, 1 extra byte for compression flag
pub const GroupG2_COMP_SIZE: usize = (2 * MODBYTES + 1) as usize;

// Byte size of element in group GT
pub const GroupGT_SIZE: usize = (12 * MODBYTES) as usize;

lazy_static! {
    pub static ref GeneratorG2: GroupG2 = GroupG2::generator();
}
