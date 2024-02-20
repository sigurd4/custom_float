
use num_traits::Signed;

use crate::{util, Dec, UInt};

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> Signed for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    U: 'static,
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:
{
    fn abs(&self) -> Self
    {
        (*self).abs()
    }

    fn abs_sub(&self, other: &Self) -> Self
    {
        (*self).abs_sub(*other)
    }

    fn signum(&self) -> Self
    {
        (*self).signum()
    }

    fn is_positive(&self) -> bool
    {
        (*self).is_sign_positive()
    }

    fn is_negative(&self) -> bool
    {
        (*self).is_sign_negative()
    }
}