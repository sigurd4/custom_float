use crate::Fp;

/// x86 extended precision floating-point format
pub type Fp80 = Fp<u128, true, 15, 1, 63, 2>;