use num_traits::Euclid;

use crate::{util, Dec, UInt};


impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> Euclid for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:
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