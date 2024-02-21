//IEEE754

#[cfg(feature = "ethnum")]
use ethnum::U256;

use crate::{util, Fp, UInt};

pub type FpHalf = Fp<u16, 5, 0, 10, 2>;
pub type FpSingle = Fp<u32, 8, 0, 23, 2>;
pub type FpDouble = Fp<u64, 11, 0, 52, 2>;
pub type FpQuadruple = Fp<u128, 15, 0, 112, 2>;
#[cfg(feature = "ethnum")]
pub type FpOctuple = Fp<U256, 19, 0, 236, 2>;

pub type DecSingle = Fp<u32, 11, 7, 13, 10>;
pub type DecDouble = Fp<u64, 13, 16, 34, 10>;
pub type DecQuadruple = Fp<u128, 17, 34, 76, 10>;