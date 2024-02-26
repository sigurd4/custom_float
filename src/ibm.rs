use crate::Fp;

/// Single precision hexadecimal floating point format
pub type HFpShort = Fp<u32, true, 7, 0, 24, 2>;

/// Double precision hexadecimal floating point format
pub type HFpLong = Fp<u64, true, 7, 0, 56, 2>;