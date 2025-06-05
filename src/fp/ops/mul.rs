use core::{num::FpCategory, ops::{Mul, MulAssign}};

use crate::fp::{UInt, Fp, util};
use super::super::as_lossless;

impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Mul<Self> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:
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
impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> MulAssign for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:
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
    use core::f32::consts::{PI, TAU};
    use std::ops::Mul;

    use test::Bencher;

    use crate::tests::{self, F};
    
    #[test]
    fn test_mul_once()
    {
        let a = F::from(PI);
        let b = F::from(TAU);
        let c = a * b;
        println!("{a} * {b} = {c}");
    }
    #[test]
    fn test_mul()
    {
        crate::tests::test_op2("mul", Mul::mul, Mul::mul, None)
    }
    #[bench]
    fn bench_mul(bencher: &mut Bencher)
    {
        test_mul();
        tests::bench_op2::<F, _>(bencher, Mul::mul)
    }
}