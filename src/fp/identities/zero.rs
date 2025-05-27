use num_traits::Zero;

use crate::{util, Fp, UInt};

impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Zero for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:
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