use core::ops::{Rem, RemAssign};

use crate::{Fp, FpRepr};
use super::super::as_lossless;

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Rem<Self> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
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
impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> RemAssign for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
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
    use core::ops::Rem;
    use test::Bencher;

    use crate::tests::{bench_op2, test_op2};

    #[test]
    fn test_rem()
    {
        test_op2!("rem", Rem::rem, None)
    }
    #[bench]
    fn bench_rem(bencher: &mut Bencher)
    {
        test_rem();
        bench_op2!(bencher, Rem::rem)
    }
}