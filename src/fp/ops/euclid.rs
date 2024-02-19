use num_traits::Euclid;

use crate::{util, Fp, UInt};


impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize> Euclid for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:
{
    #[inline]
    fn div_euclid(&self, rhs: &Self) -> Self
    {
        (*self).div_euclid(*rhs)
    }

    #[inline]
    fn rem_euclid(&self, rhs: &Self) -> Self
    {
        (*self).rem_euclid(*rhs)
    }
}