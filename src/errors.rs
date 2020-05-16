use super::types::BigNum;
use core::fmt;

#[derive(Clone, Copy)]
pub enum ValueError {
    UnequalSizeVectors(usize, usize),
    IncorrectSize(usize),
    NonPowerOf2(usize),
    OutOfRange(usize),
    NegativeValue(BigNum),
}

impl fmt::Debug for ValueError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            ValueError::UnequalSizeVectors(a, b) => ValueError::UnequalSizeVectors(*a, *b).fmt(f),
            ValueError::IncorrectSize(a) => ValueError::IncorrectSize(*a).fmt(f),
            ValueError::NonPowerOf2(a) => ValueError::NonPowerOf2(*a).fmt(f),
            ValueError::OutOfRange(a) => ValueError::OutOfRange(*a).fmt(f),
            ValueError::NegativeValue(a) => {
                write!(f, "ValueError::NegativeValue({})", a.tostring())
            }
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum SerzDeserzError {
    FieldElementBytesIncorrectSize(usize, usize),
    G1BytesIncorrectSize(usize),
    G2BytesIncorrectSize(usize),
    GTBytesIncorrectSize(usize, usize),
    RequiredHexChar,
    CannotParseFP,
    CannotParseFP2,
    CannotParseFP4,
    CannotParseG1,
    CannotParseG2,
    CannotParseGT,
}

impl fmt::Display for SerzDeserzError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            SerzDeserzError::FieldElementBytesIncorrectSize(a, b) => write!(
                f,
                "Incorrect bytes size for field element. Given {} but expected {}",
                a, b
            ),
            SerzDeserzError::G1BytesIncorrectSize(a) => write!(
                f,
                "Incorrect bytes size for G1 group element. Given {} bytes",
                a
            ),
            SerzDeserzError::G2BytesIncorrectSize(a) => write!(
                f,
                "Incorrect bytes size for G2 group element. Given {} bytes",
                a
            ),
            SerzDeserzError::GTBytesIncorrectSize(a, b) => write!(
                f,
                "Incorrect bytes size for GT group element. Given {} but expected {}",
                a, b
            ),
            SerzDeserzError::RequiredHexChar => write!(f, "Required hex character"),
            SerzDeserzError::CannotParseFP => write!(f, "Error while parsing FP"),
            SerzDeserzError::CannotParseFP2 => write!(f, "Error while parsing FP2"),
            SerzDeserzError::CannotParseFP4 => write!(f, "Error while parsing FP4"),
            SerzDeserzError::CannotParseG1 => write!(f, "Error while parsing G1"),
            SerzDeserzError::CannotParseG2 => write!(f, "Error while parsing G2"),
            SerzDeserzError::CannotParseGT => write!(f, "Error while parsing GT"),
        }
    }
}

#[macro_export]
macro_rules! check_vector_size_for_equality {
    ( $a:expr, $b:expr ) => {{
        if $a.len() != $b.len() {
            Err(ValueError::UnequalSizeVectors($a.len(), $b.len()))
        } else {
            Ok(())
        }
    }};
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_vector_size_equality() {
        let a1 = vec![1, 4, 6, 8];
        let a2 = vec![4, 5, 2, 1];
        assert!(check_vector_size_for_equality!(a1, a2).is_ok());

        let a3 = vec![1, 4, 6];
        assert!(check_vector_size_for_equality!(a3, a2).is_err());
    }
}
