use num_traits::One;

use super::{bitsize_of, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> One for Fp<U, EXP_SIZE, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 1]:
{
    fn one() -> Self
    {
        let bias = Self::exp_bias();
        Self::from_bits(bias << Self::EXP_POS)
    }
}