use crate::Fp;

/// Radeon R300 / R420 FP24
pub type Fp24 = Fp<u32, true, 7, 0, 16, 2>;