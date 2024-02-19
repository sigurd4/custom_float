//IEEE754

#[cfg(feature = "ethnum")]
use ethnum::U256;

use crate::{util, Fp, UInt};

pub type FpHalf = Fp<u16, 5, 0, 10>;
pub type FpSingle = Fp<u32, 8, 0, 23>;
pub type FpDouble = Fp<u64, 11, 0, 52>;
pub type FpQuadruple = Fp<u128, 15, 0, 112>;
#[cfg(feature = "ethnum")]
pub type FpOctuple = Fp<U256, 19, 0, 236>;

/*pub type DecSingle = Dec<u32, 11, 19>;
pub type DecDouble = Dec<u64, 13, 49>;
pub type DecQuadruple = Dec<u128, 17, 109>;*/

impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize> From<f32> for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:
{
    fn from(f: f32) -> Self
    {
        Self::from_fp::<u32, 8, 0, 23>(FpSingle::from_bits(f.to_bits()))
    }
}
impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize> Into<f32> for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:
{
    fn into(self) -> f32
    {
        f32::from_bits(FpSingle::from_fp(self).to_bits())
    }
}

impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize> From<f64> for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:
{
    fn from(f: f64) -> Self
    {
        Self::from_fp::<u64, 11, 0, 52>(FpDouble::from_bits(f.to_bits()))
    }
}
impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize> Into<f64> for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:
{
    fn into(self) -> f64
    {
        f64::from_bits(FpDouble::from_fp(self).to_bits())
    }
}

/*impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> From<f32> for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:
{
    fn from(f: f32) -> Self
    {
        Self::from_fp::<u32, 8, 0, 23>(FpSingle::from_bits(f.to_bits()))
    }
}
impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> Into<f32> for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:
{
    fn into(self) -> f32
    {
        f32::from_bits(FpSingle::from_dec(self).to_bits())
    }
}

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> From<f64> for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:
{
    fn from(f: f64) -> Self
    {
        Self::from_fp::<u64, 11, 0, 52>(FpDouble::from_bits(f.to_bits()))
    }
}
impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> Into<f64> for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:
{
    fn into(self) -> f64
    {
        f64::from_bits(FpDouble::from_dec(self).to_bits())
    }
}*/