//IEEE754

use ethnum::U256;

use crate::{bitsize_of, Fp, UInt};

pub type FpHalf = Fp<u16, 5, 10>;
pub type FpSingle = Fp<u32, 8, 23>;
pub type FpDouble = Fp<u64, 11, 52>;
pub type FpQuadruple = Fp<u128, 15, 112>;
#[cfg(feature = "ethnum")]
pub type FpOctuple = Fp<U256, 19, 236>;

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> From<f32> for Fp<U, EXP_SIZE, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 1]:
{
    fn from(f: f32) -> Self
    {
        Self::from_fp::<u32, 8, 23>(FpSingle::from_bits(f.to_bits()))
    }
}
impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> Into<f32> for Fp<U, EXP_SIZE, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 1]:
{
    fn into(self) -> f32
    {
        f32::from_bits(FpSingle::from_fp(self).to_bits())
    }
}

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> From<f64> for Fp<U, EXP_SIZE, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 1]:
{
    fn from(f: f64) -> Self
    {
        Self::from_fp::<u64, 11, 52>(FpDouble::from_bits(f.to_bits()))
    }
}
impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> Into<f64> for Fp<U, EXP_SIZE, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 1]:
{
    fn into(self) -> f64
    {
        f64::from_bits(FpDouble::from_fp(self).to_bits())
    }
}