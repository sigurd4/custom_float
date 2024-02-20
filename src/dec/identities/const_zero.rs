use num_traits::ConstZero;

use crate::{util, Dec, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> ConstZero for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:,
    U: ConstZero
{
    const ZERO: Self = Self::ZERO;
}