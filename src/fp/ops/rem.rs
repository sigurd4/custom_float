use core::ops::{Rem, RemAssign};

use crate::{UInt, Fp, util};
use super::super::as_lossless;

impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Rem<Self> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - 0 - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:
{
    type Output = Self;

    fn rem(self, rhs: Self) -> Self::Output
    {
        as_lossless!(
            [self, rhs],
            |[lhs, rhs]| [lhs % rhs],
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
                let mut divisor = Self::from_bits(Self::shift_exp(expo_a) + Self::shift_frac(frac_b));
                while divisor < dividend.divb()
                {
                    divisor = divisor.mulb();
                }
                /* compute quotient one bit at a time */
                while divisor >= fb && !divisor.is_zero() && divisor.is_finite()
                {
                    while dividend >= divisor
                    {
                        dividend -= divisor;
                    }
                    divisor = divisor.divb();
                }
                assert!(dividend.is_finite());
                /* dividend now represents remainder */
                dividend.copysign(self)
            }
        )
    }
}
impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> RemAssign for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - 0 - FRAC_SIZE]:,
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
        crate::tests::test_op2("rem", Rem::rem, Rem::rem, None)
    }
}