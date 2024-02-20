use num_traits::One;

use crate::{util, Dec, UInt};

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> One for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:
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