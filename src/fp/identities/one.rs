use num_traits::One;

use crate::{util, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> One for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:,
    [(); EXP_BASE - 2]:
{
    #[inline]
    fn one() -> Self
    {
        Self::one()
    }

    #[inline]
    fn is_one(&self) -> bool
    {
        (*self).is_one()
    }

    #[inline]
    fn set_one(&mut self)
    {
        self.set_one()
    }
}