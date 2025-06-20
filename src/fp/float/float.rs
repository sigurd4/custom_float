use core::num::FpCategory;

use num_traits::{float::FloatCore, Float};

use crate::{Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Float for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
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
        #[allow(deprecated)]
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
    use std::process::Termination;

    use num_traits::Float;
    use test::Bencher;

    use crate::{ieee754::FpDouble, tests::{bench_op1, bench_op2, test_op1, test_op2}, Fp};

    #[test]
    fn test_copysign()
    {
        test_op2!("copysign", Float::copysign, None);
    }
    #[bench]
    fn bench_copysign(bencher: &mut Bencher) -> impl Termination
    {
        test_copysign();
        bench_op2!(bencher, Fp::copysign)
    }

    #[test]
    fn test_exp2()
    {
        test_op1!("exp2", Float::exp2, None, Some(-5.0..5.0));
    }
    #[bench]
    fn bench_exp2(bencher: &mut Bencher) -> impl Termination
    {
        test_exp2();
        bench_op1!(bencher, Fp::exp2)
    }

    #[test]
    fn test_exp10()
    {
        test_op1!("exp10", |x| 10f32.powf(x), Fp::exp10, None, Some(-5.0..5.0));
    }
    #[bench]
    fn bench_exp10(bencher: &mut Bencher) -> impl Termination
    {
        test_exp10();
        bench_op1!(bencher, Fp::exp10)
    }

    #[test]
    fn test_exp()
    {
        test_op1!("exp", Float::exp, None, Some(-5.0..5.0));
    }
    #[bench]
    fn bench_exp(bencher: &mut Bencher) -> impl Termination
    {
        test_exp();
        bench_op1!(bencher, Fp::exp)
    }

    #[test]
    fn test_exp_m1()
    {
        test_op1!("exp_m1", Float::exp_m1, None, Some(-5.0..5.0));
    }
    #[bench]
    fn bench_exp_m1(bencher: &mut Bencher) -> impl Termination
    {
        test_exp_m1();
        bench_op1!(bencher, Fp::exp_m1)
    }

    #[test]
    fn test_log2()
    {
        test_op1!("log2", Float::log2, None, Some(0.5..5.0));
    }
    #[bench]
    fn bench_log2(bencher: &mut Bencher) -> impl Termination
    {
        test_log2();
        bench_op1!(bencher, Fp::log2)
    }

    #[test]
    fn test_log10()
    {
        test_op1!("log10", Float::log10, None, Some(0.5..5.0));
    }
    #[bench]
    fn bench_log10(bencher: &mut Bencher) -> impl Termination
    {
        test_log10();
        bench_op1!(bencher, Fp::log10)
    }

    #[test]
    fn test_ln()
    {
        test_op1!("ln", Float::ln, None, Some(0.5..5.0));
    }
    #[bench]
    fn bench_ln(bencher: &mut Bencher) -> impl Termination
    {
        test_ln();
        bench_op1!(bencher, Fp::ln)
    }

    #[test]
    fn test_ln_1p()
    {
        test_op1!("ln_1p", Float::ln_1p, None, Some(-0.5..5.0));
    }
    #[bench]
    fn bench_ln_1p(bencher: &mut Bencher) -> impl Termination
    {
        test_ln_1p();
        bench_op1!(bencher, Fp::ln_1p)
    }

    #[test]
    fn test_sinh()
    {
        test_op1!("sinh", Float::sinh, None, Some(-5.0..5.0));
    }
    #[bench]
    fn bench_sinh(bencher: &mut Bencher) -> impl Termination
    {
        test_sinh();
        bench_op1!(bencher, Fp::sinh)
    }

    #[test]
    fn test_cosh()
    {
        test_op1!("cosh", Float::cosh, None, Some(-5.0..5.0));
    }
    #[bench]
    fn bench_cosh(bencher: &mut Bencher) -> impl Termination
    {
        test_cosh();
        bench_op1!(bencher, Fp::cosh)
    }

    #[test]
    fn test_tanh()
    {
        test_op1!("tanh", Float::tanh, None, Some(-5.0..5.0));
    }
    #[bench]
    fn bench_tanh(bencher: &mut Bencher) -> impl Termination
    {
        test_tanh();
        bench_op1!(bencher, Fp::tanh)
    }

    #[test]
    fn test_asinh()
    {
        test_op1!("asinh", Float::asinh, None, Some(-5.0..5.0));
    }
    #[bench]
    fn bench_asinh(bencher: &mut Bencher) -> impl Termination
    {
        test_asinh();
        bench_op1!(bencher, Fp::asinh)
    }

    #[test]
    fn test_acosh()
    {
        test_op1!("acosh", Float::acosh, None, Some(1.0..5.0));
    }
    #[bench]
    fn bench_acosh(bencher: &mut Bencher) -> impl Termination
    {
        test_acosh();
        bench_op1!(bencher, Fp::acosh)
    }

    #[test]
    fn test_atanh()
    {
        test_op1!("atanh", Float::atanh, None, Some(-0.99..0.99));
    }
    #[bench]
    fn bench_atanh(bencher: &mut Bencher) -> impl Termination
    {
        test_atanh();
        bench_op1!(bencher, Fp::atanh)
    }

    #[test]
    fn test_powf_once()
    {
        type F = FpDouble;

        let x = F::from(10);
        let n = F::from(15);
        let y = x.powf(n);
        println!("{x} ** {n} = {y}");
        let xx = f32::from(x);
        let nn = f32::from(n);
        let yy = xx.powf(nn);
        if y != F::from(yy)
        {
            println!("Wrong!\n{xx} ** {nn} = {yy}");
        }
    }
    #[test]
    fn test_powf()
    {
        test_op2!("powf", Float::powf, Some(1.0));
    }
    #[bench]
    fn bench_powf(bencher: &mut Bencher) -> impl Termination
    {
        test_powf();
        bench_op2!(bencher, Fp::powf)
    }

    #[test]
    fn test_sin()
    {
        test_op1!("sin", Float::sin, None, Some(-5.0..10.0));
    }
    #[bench]
    fn bench_sin(bencher: &mut Bencher) -> impl Termination
    {
        test_sin();
        bench_op1!(bencher, Fp::sin)
    }
    
    #[test]
    fn test_cos()
    {
        test_op1!("cos", Float::cos, None, Some(-5.0..10.0));
    }
    #[bench]
    fn bench_cos(bencher: &mut Bencher) -> impl Termination
    {
        test_cos();
        bench_op1!(bencher, Fp::cos)
    }
    
    #[test]
    fn test_tan()
    {
        test_op1!("tan", Float::tan, None, Some(-1.5..10.0));
    }
    #[bench]
    fn bench_tan(bencher: &mut Bencher) -> impl Termination
    {
        test_tan();
        bench_op1!(bencher, Fp::tan)
    }
    
    #[test]
    fn test_asin()
    {
        test_op1!("asin", Float::asin, None, Some(-1.0..1.0));
    }
    #[bench]
    fn bench_asin(bencher: &mut Bencher) -> impl Termination
    {
        test_asin();
        bench_op1!(bencher, Fp::asin)
    }
    
    #[test]
    fn test_acos()
    {
        //println!("{}", FpDouble::from(0.5).acos());
        test_op1!("acos", Float::acos, None, Some(-1.0..1.0));
    }
    #[bench]
    fn bench_acos(bencher: &mut Bencher) -> impl Termination
    {
        test_acos();
        bench_op1!(bencher, Fp::acos)
    }
    
    #[test]
    fn test_atan()
    {
        test_op2!("atan", |a, b| Float::atan(a/b), Some(0.01));
        test_op1!("atan", Float::atan, Some(0.01), Some(-5.0..5.0));
    }
    #[bench]
    fn bench_atan(bencher: &mut Bencher) -> impl Termination
    {
        test_atan();
        bench_op1!(bencher, Fp::atan)
    }

    #[test]
    fn test_atan2()
    {
        test_op2!("atan2", Float::atan2, Some(0.01));
    }
    #[bench]
    fn bench_atan2(bencher: &mut Bencher) -> impl Termination
    {
        test_atan2();
        bench_op2!(bencher, Fp::atan2)
    }

    #[test]
    fn test_round()
    {
        test_op1!("round", Float::round, None, Some(-5.0..5.0));
    }
    #[bench]
    fn bench_round(bencher: &mut Bencher) -> impl Termination
    {
        test_round();
        bench_op1!(bencher, Fp::round)
    }

    #[test]
    fn test_ceil()
    {
        test_op1!("ceil", Float::ceil, None, Some(-5.0..5.0));
    }
    #[bench]
    fn bench_ceil(bencher: &mut Bencher) -> impl Termination
    {
        test_ceil();
        bench_op1!(bencher, Fp::ceil)
    }

    #[test]
    fn test_floor()
    {
        test_op1!("floor", Float::floor, None, Some(-5.0..5.0));
    }
    #[bench]
    fn bench_floor(bencher: &mut Bencher) -> impl Termination
    {
        test_floor();
        bench_op1!(bencher, Fp::floor)
    }

    #[test]
    fn test_trunc()
    {
        test_op1!("trunc", Float::trunc, None, Some(-5.0..5.0));
    }
    #[bench]
    fn bench_trunc(bencher: &mut Bencher) -> impl Termination
    {
        test_trunc();
        bench_op1!(bencher, Fp::trunc)
    }

    #[test]
    fn test_round_ties_even()
    {
        test_op1!("round_ties_even", f32::round_ties_even, Fp::round_ties_even, None, Some(-5.0..5.0));
    }
    #[bench]
    fn bench_round_ties_even(bencher: &mut Bencher) -> impl Termination
    {
        test_round_ties_even();
        bench_op1!(bencher, Fp::round_ties_even)
    }

    #[test]
    fn test_sqrt()
    {
        test_op1!("sqrt", Float::sqrt, Some(0.01), Some(0.0..5.0));
    }
    #[bench]
    fn bench_sqrt(bencher: &mut Bencher) -> impl Termination
    {
        test_sqrt();
        bench_op1!(bencher, Fp::sqrt)
    }
    
    #[test]
    fn test_cbrt()
    {
        test_op1!("cbrt", Float::cbrt, Some(0.001), Some(-5.0..5.0));
    }
    #[bench]
    fn bench_cbrt(bencher: &mut Bencher) -> impl Termination
    {
        test_cbrt();
        bench_op1!(bencher, Fp::cbrt)
    }
}