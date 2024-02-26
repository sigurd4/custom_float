use crate::Fp;

/// bfloat8 (brain float), with 8-bits using the exponent of FP16
pub type Bf8 = Fp<u8, true, 5, 0, 2, 2>;

/// bfloat16 (brain float), with 16-bits using the exponent of FP32
pub type Bf16 = Fp<u16, true, 8, 0, 7, 2>;

/// bfloat32 (brain float), with 32-bits using the exponent of FP64
pub type Bf32 = Fp<u32, true, 11, 0, 20, 2>;

/// bfloat64 (brain float), with 64-bits using the exponent of FP128
pub type Bf64 = Fp<u64, true, 15, 0, 48, 2>;

/// bfloat128 (brain float), with 128-bits using the exponent of FP256
pub type Bf128 = Fp<u128, true, 19, 0, 108, 2>;