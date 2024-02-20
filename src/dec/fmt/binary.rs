use std::fmt::Binary;

use crate::{util, Dec, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> Binary for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    U: Binary,
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:
{
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result
    {
        <U as Binary>::fmt(&self.to_bits(), f)
    }
}