use std::ops::{Rem, RemAssign};

use num_traits::ToPrimitive;

use crate::{util, Dec, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> Rem<Self> for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:
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
        let frac_b = fb.frac_bits() + (fb.int_bits() << FRAC_SIZE);
        let mut divisor = Self::from_bits((expo_a << Self::EXP_POS) + frac_b);
        if divisor.is_zero()
        {
            return Self::qnan()
        }
        while divisor < dividend.div_10()
        {
            divisor = divisor.mul_10();
        }
        /* compute quotient one bit at a time */
        while divisor >= fb
        {
            while dividend >= divisor
            {
                dividend -= divisor;
            }
            divisor = divisor.div_10();
        }
        /* dividend now represents remainder */
        dividend.copysign(self)
    }
}
impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> RemAssign for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:
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

    use crate::ieee754::DecSingle;

    #[test]
    fn test_rem()
    {
        crate::tests::test_op2_dec(Rem::rem, Rem::rem)
    }
}