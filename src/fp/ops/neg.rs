use core::ops::Neg;

use crate::{util, Fp, UInt};

impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Neg for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - 0 - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:
{
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output
    {
        if !SIGN_BIT
        {
            if self.is_zero()
            {
                return self
            }
            return Self::qnan()
        }
        Self::from_bits(self.to_bits().bitxor(Self::shift_sign(U::one())))
    }
}

#[cfg(test)]
mod test
{
    use core::ops::Neg;

    #[test]
    fn test_neg()
    {
        crate::tests::test_op1("neg", Neg::neg, Neg::neg, None, Some(-10.0..10.0))
    }
}