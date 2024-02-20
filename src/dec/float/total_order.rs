use num_traits::float::TotalOrder;

use crate::{util, Dec, UInt};

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> TotalOrder for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:
{
    #[inline]
    fn total_cmp(&self, other: &Self) -> std::cmp::Ordering
    {
        (*self).total_cmp(*other)
    }
}