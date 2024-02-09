use std::ops::{Neg, Sub, SubAssign};

use crate::{UInt, Fp, bitsize_of};

impl<U: UInt, const EXP_SIZE: usize, const INT_BIT: bool, const FRAC_SIZE: usize> Sub<Self> for Fp<U, EXP_SIZE, INT_BIT, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - INT_BIT as usize - FRAC_SIZE - 1]:,
    [(); bitsize_of::<U>() - EXP_SIZE - false as usize - FRAC_SIZE - 1]:
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output
    {
        self + (-rhs)
    }
}
impl<U: UInt, const EXP_SIZE: usize, const INT_BIT: bool, const FRAC_SIZE: usize> SubAssign for Fp<U, EXP_SIZE, INT_BIT, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - INT_BIT as usize - FRAC_SIZE - 1]:,
    [(); bitsize_of::<U>() - EXP_SIZE - false as usize - FRAC_SIZE - 1]:
{
    fn sub_assign(&mut self, rhs: Self)
    {
        *self = *self - rhs
    }
}

impl<U: UInt, const EXP_SIZE: usize, const INT_BIT: bool, const FRAC_SIZE: usize> Neg for Fp<U, EXP_SIZE, INT_BIT, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - INT_BIT as usize - FRAC_SIZE - 1]:,
    [(); bitsize_of::<U>() - EXP_SIZE - false as usize - FRAC_SIZE - 1]:
{
    type Output = Self;

    fn neg(self) -> Self::Output
    {
        Self::from_bits(self.to_bits().bitxor(U::one() << Self::BIT_SIZE - 1))
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