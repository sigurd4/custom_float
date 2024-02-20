use num_traits::{MulAdd, MulAddAssign};

use crate::{util, Dec, UInt};

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> MulAdd for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:
{
    type Output = Self;

    #[inline]
    fn mul_add(self, a: Self, b: Self) -> Self::Output
    {
        self.mul_add(a, b)
    }
}
impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> MulAddAssign for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:
{
    #[inline]
    fn mul_add_assign(&mut self, a: Self, b: Self)
    {
        *self = (*self).mul_add(a, b)
    }
}