
#[cfg(feature = "ethnum")]
use ethnum::U256;

use crate::Fp;

pub type Tf19 = Fp<u32, 8, 0, 10, 2>;
pub type Tf35 = Fp<u64, 11, 0, 23, 2>;
pub type Tf68 = Fp<u128, 15, 0, 52, 2>;
#[cfg(feature = "ethnum")]
pub type Tf132 = Fp<U256, 19, 0, 112, 2>;