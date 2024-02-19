use num_traits::Zero;

use crate::{util, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize> Zero for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:
{
    #[inline]
    fn zero() -> Self
    {
        Self::zero()
    }

    #[inline]
    fn is_zero(&self) -> bool
    {
        (*self).is_zero()
    }

    #[inline]
    fn set_zero(&mut self)
    {
        self.set_zero()
    }
}