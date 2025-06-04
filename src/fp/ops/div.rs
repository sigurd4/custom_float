use core::{num::FpCategory, ops::{Div, DivAssign}};

use crate::fp::{UInt, Fp, util};
use super::super::as_lossless;

impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Div<Self> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:
{
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output
    {
        as_lossless!(
            [self, rhs],
            |[lhs, rhs]| [lhs/rhs],
            {
                let s = self.is_sign_negative()^rhs.is_sign_negative();
                match (self.classify(), rhs.classify())
                {
                    (FpCategory::Nan, _) | (_, FpCategory::Nan) => self.add_nan(rhs).with_sign(s),
                    (FpCategory::Zero, FpCategory::Zero) | (FpCategory::Infinite, FpCategory::Infinite) => Self::qnan().with_sign(s),
                    (FpCategory::Infinite, _) | (_, FpCategory::Zero) => Self::infinity().with_sign(s),
                    (_, FpCategory::Infinite) | (FpCategory::Zero, _) => Self::zero().with_sign(s),
                    (FpCategory::Normal | FpCategory::Subnormal, FpCategory::Normal | FpCategory::Subnormal) => {
                        if rhs.abs().is_one()
                        {
                            return self.with_sign(s)
                        }
                        if self.abs().is_one()
                        {
                            return rhs.recip().with_sign(s)
                        }
                
                        let e0: U = self.exp_bits();
                        let e1: U = rhs.exp_bits();
                
                        let mut f0: U = self.mantissa_bits();
                        let mut f1: U = rhs.mantissa_bits();
                
                        if e0 == e1 && f0 == f1
                        {
                            return Self::one().with_sign(s)
                        }
                
                        let mut e = match Self::exponent_sub(e0, e1, &mut f0, &mut f1)
                        {
                            Ok(e) => e,
                            Err(done) => return done.with_sign(s)
                        };
                        let mut o = U::zero();
                        let mut f = match Self::mantissa_div(f0, f1, &mut e, &mut o)
                        {
                            Ok(f) => f,
                            Err(done) => return done.with_sign(s)
                        };
                        Self::normalize_mantissa_up(&mut e, &mut f, Some(o));
                        let mut e = match e.checked_sub(&o)
                        {
                            Some(e) => e,
                            None => return Self::zero().with_sign(s)
                        };
                
                        Self::normalize_mantissa(&mut e, &mut f, None);
                        Self::from_exp_mantissa(e, f).with_sign(s)
                    }
                }
            }
        )
    }
}
impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> DivAssign for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:
{
    #[inline]
    fn div_assign(&mut self, rhs: Self)
    {
        *self = *self / rhs
    }
}

#[cfg(test)]
mod test
{
    use std::ops::Div;

    use test::Bencher;

    use crate::tests::{self, F};

    #[test]
    fn test_div_once()
    {
        let a = F::from(4.0);
        let b = F::from(2.0);
        let c = a / b;
        println!("{a} / {b} = {c}");
    }

    #[test]
    fn test_div()
    {
        crate::tests::test_op2("div", Div::div, Div::div, None)
    }
    #[bench]
    fn bench_div(bencher: &mut Bencher)
    {
        test_div();
        tests::bench_op2::<F, _>(bencher, Div::div)
    }
}