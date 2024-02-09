use num_traits::One;

use crate::{bitsize_of, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const INT_BIT: bool, const FRAC_SIZE: usize> One for Fp<U, EXP_SIZE, INT_BIT, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - INT_BIT as usize - FRAC_SIZE - 1]:,
    [(); bitsize_of::<U>() - EXP_SIZE - false as usize - FRAC_SIZE - 1]:
{
    fn one() -> Self
    {
        let bias = Self::exp_bias();
        if INT_BIT
        {
            return Self::from_bits((bias << Self::EXP_POS) + (U::one() << Self::INT_POS))
        }
        Self::from_bits(bias << Self::EXP_POS)
    }
}