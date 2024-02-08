use num_traits::ToPrimitive;

use super::{bitsize_of, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> ToPrimitive for Fp<U, EXP_SIZE, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 1]:
{
    fn to_i64(&self) -> Option<i64>
    {
        Into::<f64>::into(*self).to_i64()
    }

    fn to_u64(&self) -> Option<u64>
    {
        Into::<f64>::into(*self).to_u64()
    }

    fn to_f32(&self) -> Option<f32>
    {
        Some(Into::<f32>::into(*self))
    }

    fn to_f64(&self) -> Option<f64>
    {
        Some(Into::<f64>::into(*self))
    }
}