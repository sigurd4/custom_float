use num_traits::float::TotalOrder;

use crate::{util, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize> PartialOrd for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering>
    {
        if self.is_nan() || other.is_nan()
        {
            return None
        }
        
        return Some(self.total_cmp(other))
    }
}