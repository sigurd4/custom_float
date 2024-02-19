use std::num::FpCategory;

use num_traits::{float::FloatCore, Float};

use crate::{util, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize> Float for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:
{
    #[inline]
    fn nan() -> Self
    {
        Self::nan()
    }

    #[inline]
    fn infinity() -> Self
    {
        Self::infinity()
    }

    #[inline]
    fn neg_infinity() -> Self
    {
        Self::neg_infinity()
    }

    #[inline]
    fn neg_zero() -> Self
    {
        Self::neg_zero()
    }

    #[inline]
    fn min_value() -> Self
    {
        Self::min_value()
    }

    #[inline]
    fn min_positive_value() -> Self
    {
        Self::min_positive_value()
    }

    #[inline]
    fn max_value() -> Self
    {
        Self::max_value()
    }

    #[inline]
    fn is_nan(self) -> bool
    {
        self.is_nan()
    }

    #[inline]
    fn is_infinite(self) -> bool
    {
        self.is_infinite()
    }

    #[inline]
    fn is_finite(self) -> bool
    {
        self.is_finite()
    }

    #[inline]
    fn is_normal(self) -> bool
    {
        self.is_normal()
    }

    #[inline]
    fn classify(self) -> FpCategory
    {
        self.classify()
    }

    #[inline]
    fn floor(self) -> Self
    {
        self.floor()
    }

    #[inline]
    fn ceil(self) -> Self
    {
        self.ceil()
    }

    #[inline]
    fn round(self) -> Self
    {
        self.round()
    }

    #[inline]
    fn trunc(self) -> Self
    {
        self.trunc()
    }

    #[inline]
    fn fract(self) -> Self
    {
        self.fract()
    }

    #[inline]
    fn abs(self) -> Self
    {
        self.abs()
    }

    #[inline]
    fn signum(self) -> Self
    {
        self.signum()
    }

    #[inline]
    fn is_sign_positive(self) -> bool
    {
        self.is_sign_positive()
    }

    #[inline]
    fn is_sign_negative(self) -> bool
    {
        self.is_sign_negative()
    }

    #[inline]
    fn mul_add(self, a: Self, b: Self) -> Self
    {
        self.mul_add(a, b)
    }

    #[inline]
    fn recip(self) -> Self
    {
        self.recip()
    }

    #[inline]
    fn powi(self, n: i32) -> Self
    {
        self.powi(n)
    }

    #[inline]
    fn powf(self, n: Self) -> Self
    {
        self.powf(n)
    }

    #[inline]
    fn sqrt(self) -> Self
    {
        self.sqrt()
    }

    #[inline]
    fn exp(self) -> Self
    {
        self.exp()
    }

    #[inline]
    fn exp2(self) -> Self
    {
        self.exp2()
    }

    #[inline]
    fn ln(self) -> Self
    {
        self.ln()
    }

    #[inline]
    fn log(self, base: Self) -> Self
    {
        self.log(base)
    }

    #[inline]
    fn log2(self) -> Self
    {
        self.log2()
    }

    #[inline]
    fn log10(self) -> Self
    {
        self.log10()
    }

    #[inline]
    fn max(self, other: Self) -> Self
    {
        self.max(other)
    }

    #[inline]
    fn min(self, other: Self) -> Self
    {
        self.min(other)
    }

    #[inline]
    fn abs_sub(self, other: Self) -> Self
    {
        self.abs_sub(other)
    }

    #[inline]
    fn cbrt(self) -> Self
    {
        self.cbrt()
    }

    #[inline]
    fn hypot(self, other: Self) -> Self
    {
        self.hypot(other)
    }

    #[inline]
    fn sin(self) -> Self
    {
        self.sin()
    }

    #[inline]
    fn cos(self) -> Self
    {
        self.cos()
    }

    #[inline]
    fn tan(self) -> Self
    {
        self.tan()
    }

    #[inline]
    fn asin(self) -> Self
    {
        self.asin()
    }

    #[inline]
    fn acos(self) -> Self
    {
        self.acos()
    }

    #[inline]
    fn atan(self) -> Self
    {
        self.atan()
    }

    #[inline]
    fn atan2(self, other: Self) -> Self
    {
        self.atan2(other)
    }

    #[inline]
    fn sin_cos(self) -> (Self, Self)
    {
        self.sin_cos()
    }

    #[inline]
    fn exp_m1(self) -> Self
    {
        self.exp_m1()
    }

    #[inline]
    fn ln_1p(self) -> Self
    {
        self.ln_1p()
    }

    #[inline]
    fn sinh(self) -> Self
    {
        self.sinh()
    }

    #[inline]
    fn cosh(self) -> Self
    {
        self.cosh()
    }

    #[inline]
    fn tanh(self) -> Self
    {
        self.tanh()
    }

    #[inline]
    fn asinh(self) -> Self
    {
        self.asinh()
    }

    #[inline]
    fn acosh(self) -> Self
    {
        self.acosh()
    }

    #[inline]
    fn atanh(self) -> Self
    {
        self.atanh()
    }

    #[inline]
    fn integer_decode(self) -> (u64, i16, i8)
    {
        <Self as FloatCore>::integer_decode(self)
    }

    #[inline]
    fn epsilon() -> Self
    {
        Self::epsilon()
    }

    #[inline]
    fn copysign(self, sign: Self) -> Self
    {
        self.copysign(sign)
    }

    #[inline]
    fn is_subnormal(self) -> bool
    {
        self.is_subnormal()
    }

    #[inline]
    fn to_degrees(self) -> Self
    {
        self.to_degrees()
    }
    
    #[inline]
    fn to_radians(self) -> Self
    {
        self.to_radians()
    }
}

#[cfg(test)]
mod test
{
    use num_traits::Float;

    use crate::{g_711::FpG711, intel::Fp80};

    #[test]
    fn test_epsilon()
    {
        let eps = FpG711::epsilon();
        assert!(eps + FpG711::one() != FpG711::one());
        assert!(eps - FpG711::min_positive_value() + FpG711::one() == FpG711::one());
    }

    #[test]
    fn test_trig()
    {
        crate::tests::test_op1(Float::acos, Float::acos);
    }

    #[test]
    fn test_exp2()
    {
        crate::tests::test_op1(Float::exp2, Float::exp2);
    }

    #[test]
    fn test_log2()
    {
        crate::tests::test_op1(Float::log2, Float::log2);
    }

    #[test]
    fn test_powf()
    {
        crate::tests::test_op2(Float::powf, Float::powf);
    }

    #[test]
    fn test_round()
    {
        crate::tests::test_op1(Float::round, Float::round);
        crate::tests::test_op1(Float::ceil, Float::ceil);
        crate::tests::test_op1(Float::floor, Float::floor);
        crate::tests::test_op1(Float::trunc, Float::trunc);
    }

    #[test]
    fn test_sqrt()
    {
        let sqrt = Fp80::from(3.0).sqrt();
        println!("{}^2 == {}", sqrt, sqrt*sqrt);
        let sqrt = (3.0).sqrt();
        println!("{}^2 == {}", sqrt, sqrt*sqrt);
    }
}