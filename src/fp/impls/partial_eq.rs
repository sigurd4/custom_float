use num_traits::Float;

use crate::{bitsize_of, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> PartialEq for Fp<U, EXP_SIZE, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 1]:
{
    fn eq(&self, other: &Self) -> bool
    {
        if self.is_nan() || other.is_nan()
        {
            return false
        }
        self.0.eq(&other.0)
    }
}