use core::cmp::Ordering;

use num_traits::float::TotalOrder;

use crate::{util, Fp, UInt};

impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> PartialOrd for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - 0 - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:
{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering>
    {
        if self.is_nan() || other.is_nan()
        {
            return None
        }
        if self.is_zero() && other.is_zero()
        {
            return Some(Ordering::Equal)
        }
        
        return Some(self.total_cmp(other))
    }
}