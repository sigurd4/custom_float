use std::ops::{Neg, Sub, SubAssign};

use crate::{util, Fp, Dec, UInt};

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> Sub<Self> for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:
{
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output
    {
        self + (-rhs)
    }
}
impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> SubAssign for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:
{
    #[inline]
    fn sub_assign(&mut self, rhs: Self)
    {
        *self = *self - rhs
    }
}

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> Neg for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:
{
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output
    {
        Self::from_bits(self.to_bits().bitxor(U::one() << Self::SIGN_POS))
    }
}

#[cfg(test)]
mod test
{
    use std::ops::Sub;

    #[test]
    fn test_sub()
    {
        crate::tests::test_op2(Sub::sub, Sub::sub)
    }
}