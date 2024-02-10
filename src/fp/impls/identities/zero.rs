use num_traits::Zero;

use crate::{bitsize_of, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const INT_BIT: bool, const FRAC_SIZE: usize> Zero for Fp<U, EXP_SIZE, INT_BIT, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - INT_BIT as usize - FRAC_SIZE - 1]:,
    [(); bitsize_of::<U>() - EXP_SIZE - false as usize - FRAC_SIZE - 1]:
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