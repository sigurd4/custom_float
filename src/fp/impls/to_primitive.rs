use num_traits::{NumCast, ToPrimitive};

use crate::{bitsize_of, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const INT_BIT: bool, const FRAC_SIZE: usize> ToPrimitive for Fp<U, EXP_SIZE, INT_BIT, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - INT_BIT as usize - FRAC_SIZE - 1]:,
    [(); bitsize_of::<U>() - EXP_SIZE - false as usize - FRAC_SIZE - 1]:
{
    fn to_i64(&self) -> Option<i64>
    {
        if self.is_sign_negative()
        {
            return (-*self).to_u64()
                .and_then(|n| <i64 as NumCast>::from(n).map(|n| -n))
        }
        self.to_u64()
            .and_then(|n| <i64 as NumCast>::from(n))
    }

    fn to_u64(&self) -> Option<u64>
    {
        self.to_uint()
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