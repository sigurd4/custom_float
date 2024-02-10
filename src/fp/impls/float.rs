use std::num::FpCategory;

use num_traits::{float::FloatCore, Float};

use crate::{bitsize_of, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const INT_BIT: bool, const FRAC_SIZE: usize> Float for Fp<U, EXP_SIZE, INT_BIT, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - INT_BIT as usize - FRAC_SIZE - 1]:,
    [(); bitsize_of::<U>() - EXP_SIZE - false as usize - FRAC_SIZE - 1]:
{
    fn nan() -> Self
    {
        Self::nan()
    }

    fn infinity() -> Self
    {
        Self::infinity()
    }

    fn neg_infinity() -> Self
    {
        Self::neg_infinity()
    }

    fn neg_zero() -> Self
    {
        Self::neg_zero()
    }

    fn min_value() -> Self
    {
        Self::min_value()
    }

    fn min_positive_value() -> Self
    {
        Self::min_positive_value()
    }

    fn max_value() -> Self
    {
        Self::max_value()
    }

    fn is_nan(self) -> bool
    {
        self.is_nan()
    }

    fn is_infinite(self) -> bool
    {
        self.is_infinite()
    }

    fn is_finite(self) -> bool
    {
        self.is_finite()
    }

    fn is_normal(self) -> bool
    {
        self.is_normal()
    }

    fn classify(self) -> FpCategory
    {
        self.classify()
    }

    fn floor(self) -> Self
    {
        self.floor()
    }

    fn ceil(self) -> Self
    {
        self.ceil()
    }

    fn round(self) -> Self
    {
        self.round()
    }

    fn trunc(self) -> Self
    {
        self.trunc()
    }

    fn fract(self) -> Self
    {
        self.fract()
    }

    fn abs(self) -> Self
    {
        self.abs()
    }

    fn signum(self) -> Self
    {
        self.signum()
    }

    fn is_sign_positive(self) -> bool
    {
        self.is_sign_positive()
    }

    fn is_sign_negative(self) -> bool
    {
        self.is_sign_negative()
    }

    fn mul_add(self, a: Self, b: Self) -> Self
    {
        self.mul_add(a, b)
    }

    fn recip(self) -> Self
    {
        self.recip()
    }

    fn powi(self, n: i32) -> Self
    {
        self.powi(n)
    }

    fn powf(self, n: Self) -> Self
    {
        self.powf(n)
    }

    fn sqrt(self) -> Self
    {
        self.sqrt()
    }

    fn exp(self) -> Self
    {
        self.exp()
    }

    fn exp2(self) -> Self
    {
        self.exp2()
    }

    fn ln(self) -> Self
    {
        self.ln()
    }

    fn log(self, base: Self) -> Self
    {
        self.log(base)
    }

    fn log2(self) -> Self
    {
        self.log2()
    }

    fn log10(self) -> Self
    {
        self.log10()
    }

    fn max(self, other: Self) -> Self
    {
        self.max(other)
    }

    fn min(self, other: Self) -> Self
    {
        self.min(other)
    }

    fn abs_sub(self, other: Self) -> Self
    {
        self.abs_sub(other)
    }

    fn cbrt(self) -> Self
    {
        self.cbrt()
    }

    fn hypot(self, other: Self) -> Self
    {
        self.hypot(other)
    }

    fn sin(self) -> Self
    {
        self.sin()
    }

    fn cos(self) -> Self
    {
        self.cos()
    }

    fn tan(self) -> Self
    {
        self.tan()
    }

    fn asin(self) -> Self
    {
        self.asin()
    }

    fn acos(self) -> Self
    {
        self.acos()
    }

    fn atan(self) -> Self
    {
        self.atan()
    }

    fn atan2(self, other: Self) -> Self
    {
        self.atan2(other)
    }

    fn sin_cos(self) -> (Self, Self)
    {
        self.sin_cos()
    }

    fn exp_m1(self) -> Self
    {
        self.exp_m1()
    }

    fn ln_1p(self) -> Self
    {
        self.ln_1p()
    }

    fn sinh(self) -> Self
    {
        self.sinh()
    }

    fn cosh(self) -> Self
    {
        self.cosh()
    }

    fn tanh(self) -> Self
    {
        self.tanh()
    }

    fn asinh(self) -> Self
    {
        self.asinh()
    }

    fn acosh(self) -> Self
    {
        self.acosh()
    }

    fn atanh(self) -> Self
    {
        self.atanh()
    }

    fn integer_decode(self) -> (u64, i16, i8)
    {
        <Self as FloatCore>::integer_decode(self)
    }

    fn epsilon() -> Self
    {
        Self::epsilon()
    }

    fn copysign(self, sign: Self) -> Self
    {
        self.copysign(sign)
    }

    fn is_subnormal(self) -> bool
    {
        self.is_subnormal()
    }

    fn to_degrees(self) -> Self
    {
        self.to_degrees()
    }
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