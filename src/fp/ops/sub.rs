use core::ops::{Sub, SubAssign};

use crate::{UInt, Fp, util};

impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Sub<Self> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:
{
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output
    {
        self.add_with_sign(rhs, true)
    }
}
impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> SubAssign for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:
{
    #[inline]
    fn sub_assign(&mut self, rhs: Self)
    {
        *self = *self - rhs
    }
}

#[cfg(test)]
mod test
{
    use std::ops::Sub;

    #[test]
    fn test_sub()
    {
        crate::tests::test_op2("sub", Sub::sub, Sub::sub, Some(0.00001))
    }
}