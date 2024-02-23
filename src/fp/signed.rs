use num_traits::Signed;

use crate::{util, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Signed for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: 'static,
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:,
    [(); EXP_BASE - 2]:
{
    fn abs(&self) -> Self
    {
        (*self).abs()
    }

    fn abs_sub(&self, other: &Self) -> Self
    {
        #[allow(deprecated)]
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