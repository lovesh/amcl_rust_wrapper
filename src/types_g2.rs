use crate::constants::MODBYTES;

// Byte size of element in group G2
pub const GroupG2_SIZE: usize = 4 * MODBYTES + 1;
pub const G2_COMP_BYTE_SIZE: usize = 2 * MODBYTES + 1;

// Byte size of element in group GT
pub const GroupGT_SIZE: usize = (12 * MODBYTES) as usize;
