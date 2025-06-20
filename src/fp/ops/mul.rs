use core::{num::FpCategory, ops::{Mul, MulAssign}};

use crate::{Fp, FpRepr};
use super::super::as_lossless;

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Mul<Self> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output
    {
        as_lossless!(
            [self, rhs],
            |[lhs, rhs]| [lhs*rhs],
            {
                // Slow! Better to avoid branching.
                /*if self.to_bits() == rhs.to_bits()
                {
                    return self.squared()
                }*/
                let s = self.is_sign_negative()^rhs.is_sign_negative();
                match (self.classify(), rhs.classify())
                {
                    (FpCategory::Nan, _) | (_, FpCategory::Nan) => self.add_nan(rhs).with_sign(s),
                    (FpCategory::Infinite, FpCategory::Zero) | (FpCategory::Zero, FpCategory::Infinite) => Self::qnan().with_sign(s),
                    (FpCategory::Zero, _) | (_, FpCategory::Zero) => Self::zero().with_sign(s),
                    (FpCategory::Infinite, _) | (_, FpCategory::Infinite) => Self::infinity().with_sign(s),
                    (FpCategory::Normal | FpCategory::Subnormal, FpCategory::Normal | FpCategory::Subnormal) => {
                        if rhs.abs().is_one()
                        {
                            return self.with_sign(s)
                        }
                        if self.abs().is_one()
                        {
                            return rhs.with_sign(s)
                        }

                        let e0: U = self.exp_bits();
                        let e1: U = rhs.exp_bits();
                
                        let mut f0: U = self.mantissa_bits();
                        let mut f1: U = rhs.mantissa_bits();
                
                        let mut e = match Self::exponent_add(e0, e1, &mut f0, Some(&mut f1))
                        {
                            Ok(e) => e,
                            Err(done) => return done.with_sign(s)
                        };
                        let mut f = match Self::mantissa_mul(f0, f1, &mut e)
                        {
                            Ok(e) => e,
                            Err(done) => return done.with_sign(s)
                        };
                
                        Self::normalize_mantissa(&mut e, &mut f, None);
                        Self::from_exp_mantissa(e, f).with_sign(s)
                    }
                }
            }
        )
    }
}
impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> MulAssign for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[inline]
    fn mul_assign(&mut self, rhs: Self)
    {
        *self = *self * rhs
    }
}

#[cfg(test)]
mod test
{
    use std::ops::Mul;

    use test::Bencher;

    use crate::{ieee754::FpHalf, tests::{bench_op2, test_op2}};
    
    #[test]
    fn test_mul_once()
    {
        type F = FpHalf;

        let a = F::one();
        let b = F::from(8388608);
        let c = a * b;
        println!("{a} * {b} = {c}");
    }
    #[test]
    fn test_mul()
    {
        test_op2!("mul", Mul::mul, Some(0.0001))
    }
    #[bench]
    fn bench_mul(bencher: &mut Bencher)
    {
        test_mul();
        bench_op2!(bencher, Mul::mul)
    }
}