use std::ops::{Rem, RemAssign};

use crate::{UInt, Fp, util};

impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Rem<Self> for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:,
    [(); EXP_BASE - 2]:
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
        let mut frac_b = fb.frac_bits();
        if !Self::IS_INT_IMPLICIT
        {
            frac_b = frac_b + (fb.int_bits() << Self::INT_POS)
        }
        let mut divisor = Self::from_bits((expo_a << Self::EXP_POS) + frac_b);
        while divisor < dividend
        {
            divisor = divisor.mul_base();
        }
        /* compute quotient one bit at a time */
        while divisor >= fb && !divisor.is_zero()
        {
            while dividend >= divisor
            {
                dividend -= divisor;
            }
            divisor = divisor.div_base();
        }
        /* dividend now represents remainder */
        dividend.copysign(self)
    }
}
impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> RemAssign for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:,
    [(); EXP_BASE - 2]:
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