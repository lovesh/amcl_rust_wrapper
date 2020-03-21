use std::ops::BitXor;

#[derive(Debug, Eq, PartialEq)]
pub enum Sgn0 {
    /// Either 0 or positive
    NonNegative,
    /// Neither 0 nor positive
    Negative
}

impl BitXor for Sgn0 {
    type Output = Self;

    fn bitxor(self, rhs: Self) -> Self {
        if self == rhs {
            Sgn0::NonNegative
        } else {
            Sgn0::Negative
        }
    }
}

