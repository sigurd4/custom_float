
#[cfg(feature = "ethnum")]
use ethnum::U256;

use crate::Fp;

/// 19-bit tensor float using the exponent of FP32 and the mantissa of FP16
pub type Tf19 = Fp<u32, true, 8, 0, 10, 2>;

/// 35-bit tensor float using the exponent of FP64 and the mantissa of FP32
pub type Tf35 = Fp<u64, true, 11, 0, 23, 2>;

/// 68-bit tensor float using the exponent of FP128 and the mantissa of FP64
pub type Tf68 = Fp<u128, true, 15, 0, 52, 2>;

/// 132-bit tensor float using the exponent of FP256 and the mantissa of FP128
#[cfg(feature = "ethnum")]
pub type Tf132 = Fp<U256, true, 19, 0, 112, 2>;