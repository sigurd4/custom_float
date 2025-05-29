#![cfg_attr(not(test), no_std)]
#![feature(test)]
#![allow(incomplete_features)]
#![feature(const_trait_impl)]
#![feature(generic_const_exprs)]
#![feature(trait_alias)]
#![feature(iter_next_chunk)]
#![feature(portable_simd)]
#![feature(float_gamma)]
#![feature(let_chains)]
#![feature(specialization)]
#![feature(f16)]
#![feature(f128)]
#![feature(generic_arg_infer)]
#![feature(assert_matches)]
#![feature(decl_macro)]
#![feature(more_float_constants)]
#![feature(duration_millis_float)]
#![allow(clippy::excessive_precision)]

//! # Custom Float
//!
//! This crate adds a custom floating point number type, `Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>`, where the bit size of the exponent and mantissa can
//! be set separately.
//!
//! `U` is the underlying unsigned integer type which is used to represent the number.
//!
//! `SIGN_BIT` is wether or not the number has a sign bit.
//!
//! `EXP_SIZE` is the size of the exponent in bits.
//!
//! `INT_SIZE` is the size of the integer part of the mantissa in bits. If zero, then the integer bit is implicit.
//!
//! `FRAC_SIZE` is the size of the fractional part of the mantissa in bits.
//!
//! `EXP_BASE` is the base of the exponent.
//!
//! The total bit size of `U` must be greater or equal to `SIGN_BIT` + `EXP_SIZE` + `INT_SIZE` + `FRAC_SIZE` to contain the entire number.
//!
//! The bit layout is as follows:
//! ```txt
//! No data: | Sign:      | Exponent:  | Integer:   | Fractional: |
//! <  ..  > | <SIGN_BIT> | <EXP_SIZE> | <INT_SIZE> | <FRAC_SIZE> |
//! ```
//!
//! The value of a real floating-point number is the following:
//! ```txt
//! x = (-1)**sign*EXP_BASE**(exponent - bias)*mantissa
//! ```
//!
//! where the bias equals
//! ```txt
//! bias = 2**(EXP_SIZE - 1) - 1
//! ```
//!
//! If the exponent has the maximum value, the number is either infinity or NaN.
//!
//! The number then automatically implements `num::Float`, and supports all ordinary floating point operations.
//!
//! This allows simple implementation of special floating point types, such as TensorFloat, IEEE754 Quadruple/binary128, Fp80, and BFloat16.
//!
//! The accuracy of all of the floating point operations are not perfect, but work well enough to be usable. Various plots showing the accuracy of basic functions are shown in the [plots](https://github.com/sigurd4/custom_float/tree/master/plots) subfolder.
//!
//! All floats can be converted into each other painlessly, though the conversion may produce rounding errors or unbounded outputs when converting to a float with lesser
//! resolution.
//!
//! # Examples
//!
//! ```rust
//! #![feature(generic_const_exprs)]
//!
//! use custom_float::Fp;
//!
//! type FpSingle = Fp<u32, true, 8, 0, 23, 2>;
//!
//! let two = FpSingle::from(2);
//! let four = FpSingle::from(4);
//!
//! assert_eq!(two + two, four);
//! ```

use num_traits::{Bounded, CheckedNeg, CheckedShl, CheckedShr, PrimInt, Signed, Unsigned};

moddef::moddef!(
    pub mod {
        ati,
        g_711,
        google,
        ibm,
        ieee754,
        intel,
        khronos,
        nvidia,
        pixar
    },
    flat(pub) mod {
        fp
    },
    mod {
        util,
        plot for cfg(test)
    }
);

#[cfg(test)]
extern crate test;

pub trait AnyInt = Bounded + PrimInt + CheckedShl + CheckedShr;
pub trait UInt = Unsigned + AnyInt + core::fmt::Debug + core::fmt::Binary;
pub trait Int = Signed + AnyInt + CheckedNeg;

#[cfg(test)]
mod tests
{
    #![allow(unused)]

    use core::{f64::consts::*, ops::Neg};
    use std::{
        ops::{Range, RangeBounds}, process::Termination, time::{Instant, SystemTime}
    };

    use linspace::LinspaceArray;
    use num::Complex;
    use num_traits::{Float, Inv, One, ToPrimitive, Zero};
    use test::Bencher;

    use crate::{
        ati::Fp24,
        g_711::FpG711,
        google::{Bf16, Bf32, Bf8},
        ibm::{HFpLong, HFpShort},
        ieee754::{DecDouble, FpDouble, FpHalf, FpQuadruple, FpSingle},
        intel::Fp80,
        khronos::KhronosFp11,
        nvidia::Tf19,
        plot, Fp
    };

    pub type F = FpDouble;

    #[test]
    fn it_works()
    {
        let mut x = [0.0; 8].map(|x| F::from(x));
        x[0] = F::one() + F::one() + F::one();
        x[1] = F::one() + F::from(0.5);
        x[2] = F::one();

        println!("{:?}", x);

        //let mut y = [Complex::zero(); 5];
        //x.real_fft(&mut y);
        //x.real_ifft(&y);

        //println!("{:?}", x);
    }

    #[test]
    fn test_dec()
    {
        let x = 1e-1;
        let d = DecDouble::from(x);

        println!("{:e}", d);

        let e = d.exp_bits() as i64 - DecDouble::exp_bias() as i64;
        let i = d.int_bits();
        let f = d.frac_bits();

        let y = (i as f64 + f as f64/2.0f64.powi(DecDouble::FRAC_SIZE as i32))*10.0f64.powi(e as i32);
        println!("{:e}", y);

        let x = x.into();

        assert!(FpDouble::from(y).approx_eq(x));

        let y = FpDouble::from_fp(d);

        assert!(y.approx_eq(x));
    }

    pub fn ttable<F: Float>() -> Vec<F>
    {
        let a = [
            F::one(),
            F::from(10000.0).unwrap(),
            F::from(0.5).unwrap(),
            F::from(3.333333).unwrap(),
            F::from(10.0).unwrap(),
            F::from(16.0).unwrap(),
            F::from(2.2).unwrap(),
            F::from(0.3643634).unwrap(),
            F::from(0.1353856).unwrap(),
            F::from(0.035136474).unwrap(),
            F::from(0.74378643).unwrap(),
            F::from(1.253464).unwrap(),
            F::from(PI).unwrap(),
            F::from(TAU).unwrap(),
            F::from(PHI).unwrap(),
            F::from(EGAMMA).unwrap(),
            F::from(FRAC_PI_2).unwrap(),
            F::from(FRAC_PI_3).unwrap(),
            F::from(FRAC_PI_4).unwrap(),
            F::from(FRAC_PI_6).unwrap(),
            F::from(FRAC_PI_8).unwrap(),
            F::from(FRAC_1_PI).unwrap(),
            F::from(FRAC_1_SQRT_PI).unwrap(),
            F::from(FRAC_1_SQRT_2PI).unwrap(),
            F::from(FRAC_2_PI).unwrap(),
            F::from(FRAC_2_SQRT_PI).unwrap(),
            F::from(SQRT_2).unwrap(),
            F::from(FRAC_1_SQRT_2).unwrap(),
            F::from(SQRT_3).unwrap(),
            F::from(FRAC_1_SQRT_3).unwrap(),
            F::from(E).unwrap(),
            F::from(LOG2_10).unwrap(),
            F::from(LOG2_E).unwrap(),
            F::from(LOG10_2).unwrap(),
            F::from(LOG10_E).unwrap(),
            F::from(LN_2).unwrap(),
            F::from(LN_10).unwrap(),
            F::epsilon(),
            F::from(1.0.tan()).unwrap(),
            F::one() + F::epsilon(),
        ];

        let mut a = a.into_iter()
            .chain(a.into_iter()
                .map(|x| x.recip())
            ).chain([
                F::min_positive_value(),
                F::zero(),
                F::infinity(),
            ]).collect::<Vec<_>>();

        a.clone()
            .into_iter()
            .chain(a.into_iter()
                .map(Neg::neg)
            ).chain([
                F::nan()
            ])
            .collect()
    }

    pub(crate) fn matches<F>(a: F, b: F, tol: Option<F>) -> bool
    where
        F: Float
    {
        (match tol
        {
            Some(tol) => !((a - b).abs() > tol),
            None => !(a != b)
        }) || a.is_nan() && b.is_nan()
    }

    pub fn test_op2(fn_name: &str, op1: impl Fn(f32, f32) -> f32, op2: impl Fn(F, F) -> F, d: Option<f32>)
    {
        for f0 in crate::tests::ttable()
        {
            for f1 in crate::tests::ttable()
            {
                let fp0 = F::from(f0);
                let fp1 = F::from(f1);

                //println!("{} ? {}", f0, f1);

                let s = op1(f0, f1);
                let sp: f32 = op2(fp0, fp1).into();

                if !matches(s, sp, d)
                {
                    if f0.is_subnormal()
                    {
                        println!("f is subnormal");
                    }
                    if s.is_subnormal()
                    {
                        println!("y is subnormal");
                    }
                    println!("{:?} ? {:?} == {:?} != {:?}", f0, f1, s, sp);
                }
            }
        }
    }

    pub fn test_op1(fn_name: &str, op1: impl Fn(f32) -> f32, op2: impl Fn(F) -> F, d: Option<f32>, r: Option<Range<f32>>)
    {
        for f0 in ttable::<f32>()
        {
            let fp0 = F::from(f0);

            let s = op1(f0);
            let sp: f32 = op2(fp0).into();

            if !matches(s, sp, d)
            {
                if f0.is_subnormal()
                {
                    println!("f is subnormal");
                }
                if s.is_subnormal()
                {
                    println!("y is subnormal");
                }
                println!("{:e} ? == {:e} != {:e}", f0, s, sp);
            }
        }

        if let Some(r) = r
        {
            #[cfg(debug_assertions)]
            plot_approx(fn_name, r.clone(), &op1, |x| op2(Fp::from(x)).into());
            #[cfg(debug_assertions)]
            plot_err(fn_name, r.clone(), &op1, |x| op2(Fp::from(x)).into());
            #[cfg(not(debug_assertions))]
            {
                const NANOS_PER_MILLIS: f64 = 1e6;
                plot_bench(
                    fn_name,
                    r,
                    |x| {
                        let t0 = Instant::now();

                        for _ in 0..M
                        {
                            let _ = op1(x);
                        }

                        Instant::now().duration_since(t0).div_f64(M as f64/NANOS_PER_MILLIS).as_millis_f32()
                    },
                    |x| {
                        let x = Fp::from(x);
                        let t0 = Instant::now();

                        for _ in 0..M
                        {
                            let _ = op2(x);
                        }

                        Instant::now().duration_since(t0).div_f64(M as f64/NANOS_PER_MILLIS).as_millis_f32()
                    }
                )
            }
        }
    }

    pub fn bench_op1<F, O>(bencher: &mut Bencher, mut op: impl FnMut(F) -> O)
    where
        F: Float
    {
        let mut x = ttable::<F>()
            .into_iter()
            .cycle();

        bencher.iter(|| op(x.next().unwrap()));
    }

    pub fn bench_op2<F, O>(bencher: &mut Bencher, mut op: impl FnMut(F, F) -> O)
    where
        F: Float
    {
        let mut x = ttable::<F>()
            .into_iter()
            .flat_map(|lhs| ttable::<F>()
                .into_iter()
                .map(move |rhs| (lhs, rhs))
            ).collect::<Vec<_>>()
            .into_iter()
            .cycle();

        bencher.iter(|| {
            let (lhs, rhs) = x.next().unwrap();
            op(lhs, rhs)
        });
    }

    const M: usize = 64;
    const N: usize = 1024;
    const PLOT_TARGET: &str = "plots";

    #[allow(unused)]
    pub fn plot_err<R>(fn_name: &str, range: R, func: impl Fn(f32) -> f32, approx: impl Fn(f32) -> f32)
    where
        R: RangeBounds<f32> + LinspaceArray<f32, N>
    {
        let x: [f32; N] = range.linspace_array();
        let y_approx = x.map(approx);

        let y = x.map(func);

        let e = y_approx.into_iter().zip(y).map(|(y_approx, y)| y_approx - y).next_chunk().unwrap();

        let plot_title: &str = &format!("{fn_name}(x) error");
        let plot_path: &str = &format!("{PLOT_TARGET}/error/{fn_name}_error.png");

        plot::plot_curves(plot_title, plot_path, [x], [e]).expect("Plot error");

        /*let (avg_error, max_abs_error) = y.zip(y_approx)
            .map(|(y, y_approx)| y - y_approx)
            .map(|y| (y, y.abs()))
            .reduce(|a, b| (a.0 + b.0, a.1.max(b.1)))
            .map(|(sum_error, max_abs_error)| (sum_error/N as f32, max_abs_error))
            .unwrap_or_default();
        println!("Average Error: {}", avg_error);
        println!("Max |Error|: {}", max_abs_error);*/
    }

    #[allow(unused)]
    pub fn plot_approx<R>(fn_name: &str, range: R, func: impl Fn(f32) -> f32, approx: impl Fn(f32) -> f32)
    where
        R: RangeBounds<f32> + LinspaceArray<f32, N>
    {
        let x: [f32; N] = range.linspace_array();
        let y_approx = x.map(approx);

        let y = x.map(func);

        let plot_title: &str = &format!("{fn_name}(x)");
        let plot_path: &str = &format!("{PLOT_TARGET}/{fn_name}.png");

        plot::plot_curves(plot_title, plot_path, [x, x], [y, y_approx]).expect("Plot error");

        /*let (avg_error, max_abs_error) = y.zip(y_approx)
            .map(|(y, y_approx)| y - y_approx)
            .map(|y| (y, y.abs()))
            .reduce(|a, b| (a.0 + b.0, a.1.max(b.1)))
            .map(|(sum_error, max_abs_error)| (sum_error/N as f32, max_abs_error))
            .unwrap_or_default();
        println!("Average Error: {}", avg_error);
        println!("Max |Error|: {}", max_abs_error);*/
    }

    #[allow(unused)]
    pub fn plot_bench<R>(fn_name: &str, range: R, func: impl Fn(f32) -> f32, approx: impl Fn(f32) -> f32)
    where
        R: RangeBounds<f32> + LinspaceArray<f32, N>
    {
        let x: [f32; N] = range.linspace_array();
        let y_approx = x.map(approx);

        let y = x.map(func);

        let plot_title: &str = &format!("{fn_name}(x) benchmark");
        let plot_path: &str = &format!("{PLOT_TARGET}/bench/{fn_name}_bench.png");

        plot::plot_curves(plot_title, plot_path, [x, x], [y, y_approx]).expect("Plot error");

        /*let (avg_error, max_abs_error) = y.zip(y_approx)
            .map(|(y, y_approx)| y - y_approx)
            .map(|y| (y, y.abs()))
            .reduce(|a, b| (a.0 + b.0, a.1.max(b.1)))
            .map(|(sum_error, max_abs_error)| (sum_error/N as f32, max_abs_error))
            .unwrap_or_default();
        println!("Average Error: {}", avg_error);
        println!("Max |Error|: {}", max_abs_error);*/
    }

    #[test]
    fn int_bit()
    {
        let x = Fp80::from(2.0);
        println!("{}", x.int_bits());
        println!("{}", x);
        println!("{:b}", x)
    }

    #[test]
    fn test_ident()
    {
        test_op1("ident", |x| x, |x| x, None, Some(-5.0..20.0))
    }
    #[bench]
    fn bench_ident(bencher: &mut Bencher) -> impl Termination
    {
        test_ident();
        bench_op1::<F, _>(bencher, |x| x);
    }

    #[test]
    fn test_to_f16()
    {
        test_op1("to_f16", |x| (x as f16) as _, |x| (f16::from(x) as f32).into(), None, Some(-5.0..20.0))
    }
    #[bench]
    fn bench_to_f16(bencher: &mut Bencher) -> impl Termination
    {
        test_to_f16();
        bench_op1::<F, _>(bencher, f16::from);
    }

    #[test]
    fn test_to_f32()
    {
        test_op1("to_f32", |x| (x as f32) as _, |x| f32::from(x).into(), None, Some(-5.0..20.0))
    }
    #[bench]
    fn bench_to_f32(bencher: &mut Bencher) -> impl Termination
    {
        test_to_f32();
        bench_op1::<F, _>(bencher, f32::from);
    }

    #[test]
    fn test_to_f64()
    {
        test_op1("to_f64", |x| (x as f64) as _, |x| (f64::from(x) as f32).into(), None, Some(-5.0..20.0))
    }
    #[bench]
    fn bench_to_f64(bencher: &mut Bencher) -> impl Termination
    {
        test_to_f64();
        bench_op1::<F, _>(bencher, f64::from);
    }

    #[test]
    fn test_to_f128()
    {
        test_op1("to_f128", |x| (x as f128) as f32, |x| (f128::from(x) as f32).into(), None, Some(-5.0..20.0))
    }
    #[bench]
    fn bench_to_f128(bencher: &mut Bencher) -> impl Termination
    {
        test_to_f128();
        bench_op1::<F, _>(bencher, f128::from);
    }

    #[test]
    fn test_to_int()
    {
        for n in u8::MIN..=u8::MAX
        {
            let f = F::from_uint(n);
            assert_eq!(f.to_uint_wrapping::<u8>(), n);
            assert_eq!(f.to_uint(), Some(n));
        }
    }

    #[test]
    fn test_show_all()
    {
        type F = Fp<u8, true, 3, 0, 4, 16>;

        let mut x: F = F::from(0);

        let mut x_prev = x;
        while x.is_finite()
        {
            println!("{}", x);
            x = x.next_up();
            if x.is_finite()
            {
                if !x.is_zero()
                {
                    assert_ne!(x, x_prev);
                }
                assert!(x > x_prev);
                assert_eq!(x_prev, x.next_down());
                x_prev = x;
            }
        }
    }

    #[test]
    fn test_next_up_down_unsigned()
    {
        type F = Fp<u8, false, 3, 1, 3, { usize::MAX }>;

        let mut x = F::zero();

        loop
        {
            let y = x.next_up();
            if !y.is_zero()
            {
                if x == y
                {
                    let y = x.next_up();
                }
                assert_ne!(x, y);
            }
            if !(x == y.next_down())
            {
                let yy = x.next_up();
                let xx = y.next_down();
                println!("{x} ^ {yy} v {xx}")
            }
            assert_eq!(x, y.next_down());
            x = y;
            if !x.is_finite()
            {
                break;
            }
        }

        loop
        {
            let y = x.next_down();
            if !y.is_zero()
            {
                assert_ne!(x, y);
            }
            assert_eq!(x, y.next_up());
            x = y;
            if !x.is_zero()
            {
                break;
            }
        }
    }
}
