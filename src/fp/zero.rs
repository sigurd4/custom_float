use num_traits::Zero;

use super::{bitsize_of, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> Zero for Fp<U, EXP_SIZE, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 1]:
{
    fn zero() -> Self
    {
        Self::from_bits(U::zero())
    }

    fn is_zero(&self) -> bool
    {
        (self.to_bits() << Self::BIT_SIZE - Self::SIGN_POS).is_zero()
    }
}