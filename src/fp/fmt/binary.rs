use core::fmt::Binary;

use crate::{util, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Binary for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: Binary,
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:,
    [(); EXP_BASE - 2]:
{
    #[inline]
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result
    {
        <U as Binary>::fmt(&self.to_bits(), f)
    }
}