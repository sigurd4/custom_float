use num_traits::One;

use crate::{bitsize_of, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const INT_BIT: bool, const FRAC_SIZE: usize> One for Fp<U, EXP_SIZE, INT_BIT, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - INT_BIT as usize - FRAC_SIZE - 1]:,
    [(); bitsize_of::<U>() - EXP_SIZE - false as usize - FRAC_SIZE - 1]:
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