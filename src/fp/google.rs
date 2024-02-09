use crate::Fp;

pub type Bf8 = Fp<u8, 5, false, 2>;
pub type Bf16 = Fp<u16, 8, false, 7>;
pub type Bf32 = Fp<u32, 11, false, 20>;
pub type Bf64 = Fp<u64, 15, false, 48>;
pub type Bf128 = Fp<u128, 19, false, 108>;