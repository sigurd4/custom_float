use crate::{ieee754::{FpDouble, FpSingle}, util, Fp, UInt};


impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Into<f32> for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:,
    [(); EXP_BASE - 2]:
{
    fn into(self) -> f32
    {
        f32::from_bits(FpSingle::from_fp(self).to_bits())
    }
}
impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Into<f64> for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:,
    [(); EXP_BASE - 2]:
{
    fn into(self) -> f64
    {
        f64::from_bits(FpDouble::from_fp(self).to_bits())
    }
}