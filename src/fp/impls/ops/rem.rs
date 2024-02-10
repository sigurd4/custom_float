use std::ops::{Rem, RemAssign};

use crate::{UInt, Fp, bitsize_of};

impl<U: UInt, const EXP_SIZE: usize, const INT_BIT: bool, const FRAC_SIZE: usize> Rem<Self> for Fp<U, EXP_SIZE, INT_BIT, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - INT_BIT as usize - FRAC_SIZE - 1]:,
    [(); bitsize_of::<U>() - EXP_SIZE - false as usize - FRAC_SIZE - 1]:
{
    type Output = Self;

    fn rem(self, rhs: Self) -> Self::Output
    {
        if self.is_nan() || rhs.is_nan()
        {
            return self + rhs;
        }
        if self.is_infinite() || rhs.is_zero()
        {
            return Self::snan()
        }

        let fa = self.abs();
        let fb = rhs.abs();

        if fa < fb
        {
            return self
        }
    
        let mut dividend = fa;
        /* normalize divisor */
        let expo_a = fa.exp_bits();
        let frac_b = fb.frac_bits();
        let mut divisor = Self::from_bits((expo_a << Self::EXP_POS) + frac_b);
        if divisor <= dividend*Self::from(0.5)
        {
            divisor += divisor;
        }
        /* compute quotient one bit at a time */
        while divisor >= fb
        {
            if dividend >= divisor
            {
                dividend -= divisor;
            }
            divisor *= Self::from(0.5);
        }
        while dividend >= fb
        {
            dividend -= fb;
        }
        /* dividend now represents remainder */
        dividend.copysign(self)
    }
}
impl<U: UInt, const EXP_SIZE: usize, const INT_BIT: bool, const FRAC_SIZE: usize> RemAssign for Fp<U, EXP_SIZE, INT_BIT, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - INT_BIT as usize - FRAC_SIZE - 1]:,
    [(); bitsize_of::<U>() - EXP_SIZE - false as usize - FRAC_SIZE - 1]:
{
    #[inline]
    fn rem_assign(&mut self, rhs: Self)
    {
        *self = *self % rhs
    }
}

#[cfg(test)]
mod test
{
    use std::ops::Rem;

    #[test]
    fn test_rem()
    {
        crate::tests::test_op2(Rem::rem, Rem::rem)
    }
}