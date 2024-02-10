use num_traits::float::TotalOrder;

use crate::{bitsize_of, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const INT_BIT: bool, const FRAC_SIZE: usize> TotalOrder for Fp<U, EXP_SIZE, INT_BIT, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - INT_BIT as usize - FRAC_SIZE - 1]:,
    [(); bitsize_of::<U>() - EXP_SIZE - false as usize - FRAC_SIZE - 1]:
{
    fn total_cmp(&self, other: &Self) -> std::cmp::Ordering
    {
        let s1 = self.sign_bit();
        let s2 = self.sign_bit();

        if s1 != s2
        {
            return s2.cmp(&s1)
        }

        let s = !s1.is_zero();
        
        let mut left = self.to_bits();
        let mut right = other.to_bits();

        left = left ^ (((left >> (bitsize_of::<U>() - 1))) >> 1usize);
        right = right ^ (((right >> (bitsize_of::<U>() - 1))) >> 1usize);

        if s {right.cmp(&left)} else {left.cmp(&right)}
    }
}