use core::num::FpCategory;

use num_traits::{MulAdd, MulAddAssign};

use crate::{util, Fp, UInt};
use super::super::as_lossless;

impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:
{
    /// Fused multiply-add. Computes `(self * a) + b`.
    ///
    /// # Examples
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let m = FpDouble::from(10.0);
    /// let x = FpDouble::from(4.0);
    /// let b = FpDouble::from(60.0);
    ///
    /// // 100.0
    /// let abs_difference = (m.mul_add(x, b) - (m*x + b)).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn mul_add(self, a: Self, b: Self) -> Self
    {
        as_lossless!(
            [self, a, b],
            |[a, b, c]| [a.mul_add(b, c)],
            {
                let s0 = self.is_sign_negative()^a.is_sign_negative();
                let s1 = b.is_sign_negative();
                match (self.classify(), a.classify())
                {
                    (FpCategory::Nan, _) | (_, FpCategory::Nan) => self.add_nan(a).with_sign(s0).add_nan(b),
                    (FpCategory::Infinite, FpCategory::Zero) | (FpCategory::Zero, FpCategory::Infinite) => Self::qnan().with_sign(s0).add_nan(b),
                    (FpCategory::Zero, _) | (_, FpCategory::Zero) => b,
                    (FpCategory::Infinite, _) | (_, FpCategory::Infinite) => match b.classify()
                    {
                        FpCategory::Nan => b,
                        FpCategory::Infinite if s0 != s1 => Self::qnan(),
                        FpCategory::Infinite | FpCategory::Normal | FpCategory::Subnormal | FpCategory::Zero => Self::infinity().with_sign(s0)
                    },
                    (FpCategory::Normal | FpCategory::Subnormal, FpCategory::Normal | FpCategory::Subnormal) => match b.classify()
                    {
                        FpCategory::Nan | FpCategory::Infinite => b,
                        FpCategory::Zero => self*a,
                        FpCategory::Normal | FpCategory::Subnormal => {
                            if a.abs().is_one()
                            {
                                return self.with_sign(s0) + b
                            }
                            if self.abs().is_one()
                            {
                                return a.with_sign(s0) + b
                            }
            
                            let e0 = self.exp_bits();
                            let e1 = a.exp_bits();
                            let e2 = b.exp_bits();
                    
                            let mut f0 = self.mantissa_bits();
                            let mut f1 = a.mantissa_bits();
                            let mut f2 = b.mantissa_bits();
                    
                            let mut e = match Self::exponent_add(e0, e1, &mut f0, Some(&mut f1))
                            {
                                Ok(e) => e,
                                Err(done) => return done.with_sign(s0) + b
                            };
                            let mut f = match Self::mantissa_mul(f0, f1, &mut e)
                            {
                                Ok(e) => e,
                                Err(done) => return done.with_sign(s0) + b
                            };

                            e = Self::max_exponents(e, e2, &mut f, &mut f2);
                            let s = Self::add_signs(s0, s1, f, f2);
                            let mut f = Self::abs_add_mantissas(&mut e, f, f2, s0 != s1);
                    
                            Self::normalize_mantissa(&mut e, &mut f, None);
                            Self::from_exp_mantissa(e, f).with_sign(s)
                        }
                    }
                }
            }
        )
    }
}

impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> MulAdd for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:
{
    type Output = Self;

    #[inline]
    fn mul_add(self, a: Self, b: Self) -> Self::Output
    {
        self.mul_add(a, b)
    }
}
impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> MulAddAssign for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:
{
    #[inline]
    fn mul_add_assign(&mut self, a: Self, b: Self)
    {
        *self = (*self).mul_add(a, b)
    }
}

#[cfg(test)]
mod test
{
    use num_traits::MulAdd;
    use test::Bencher;

    use crate::tests::{self, F};
    
    #[test]
    fn test_mul_add_once()
    {
        let a = F::from(2f32);
        let b = F::from(2f32);
        let c = F::from(2f32);
        let d = a.mul_add(b, c);
        println!("{a} * {b} + {c} = {d}");
    }
    #[test]
    fn test_mul_add()
    {
        crate::tests::test_op3("mul_add", MulAdd::mul_add, MulAdd::mul_add, Some(0.001))
    }
    #[bench]
    fn bench_mul_add(bencher: &mut Bencher)
    {
        test_mul_add();
        tests::bench_op3::<F, _>(bencher, F::mul_add)
    }
}