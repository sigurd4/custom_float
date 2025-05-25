#![cfg_attr(not(test), no_std)]
#![allow(incomplete_features)]
#![feature(const_trait_impl)]
#![feature(generic_const_exprs)]
#![feature(trait_alias)]
#![feature(iter_next_chunk)]
#![feature(portable_simd)]
#![feature(float_gamma)]

//! # Custom Float
//! 
//! This crate adds a custom floating point number type, `Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>`, where the bit size of the exponent and mantissa can be set separately.
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
//! All floats can be converted into each other painlessly, though the conversion may produce rounding errors or unbounded outputs when converting to a float with lesser resolution.
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

use core::fmt::Debug;

use num_traits::{CheckedNeg, CheckedShl, CheckedShr, PrimInt, Signed, Unsigned};

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

pub trait UInt = Unsigned + PrimInt + CheckedShl + CheckedShr + Debug;
pub trait Int = Signed + PrimInt + CheckedShl + CheckedShr + CheckedNeg;

#[cfg(test)]
mod tests
{
    #![allow(unused)]

    use std::{ops::{Range, RangeBounds}, time::{Instant, SystemTime}};

    use linspace::LinspaceArray;
    use num::Complex;
    use num_traits::{Float, One, ToPrimitive, Zero};

    use crate::{ati::Fp24, g_711::FpG711, google::{Bf16, Bf32, Bf8}, ibm::{HFpLong, HFpShort}, ieee754::{DecDouble, FpDouble, FpHalf, FpQuadruple, FpSingle}, intel::Fp80, khronos::KhronosFp11, nvidia::Tf19, plot, Fp};

    pub type F = Fp<u32, true, 8, 0, 23, 2>;

    #[test]
    fn it_works()
    {
        let mut x = [0.0; 8].map(|x| F::from(x));
        x[0] = F::one() + F::one() + F::one();
        x[1] = F::one() + F::one();
        x[2] = F::one();

        println!("{:?}", x);

        //let mut y = [Complex::zero(); 5];
        //x.real_fft(&mut y);
        //x.real_ifft(&y);

        //println!("{:?}", x);
    }

    pub fn ttable<F: Float>() -> Vec<F>
    {
        vec![
            F::from(3.333333).unwrap(),
            F::from(10.0).unwrap(),
            F::from(16.0).unwrap(),
            F::from(-2.2).unwrap(),
            F::from(2.2).unwrap(),
            F::one(),
            F::one().tan(),
            F::one() + F::epsilon(),
            F::zero(),
            F::from(0.5).unwrap(),
            F::from(-0.5).unwrap(),
            F::nan(),
            F::infinity(),
            F::neg_infinity(),
            F::min_positive_value(),
        ]
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

                if match d
                {
                    Some(d) => (s - sp).abs() > d,
                    None => s != sp
                } && !(s.is_nan() && sp.is_nan())
                {
                    if f0.is_subnormal()
                    {
                        println!("f is subnormal");
                    }
                    if s.is_subnormal()
                    {
                        println!("y is subnormal");
                    }
                    println!("{:e} ? {:e} == {:e} != {:e}", f0, f1, s, sp);
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

            if match d
            {
                Some(d) => (s - sp).abs() > d,
                None => s != sp
            } && !(s.is_nan() && sp.is_nan())
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
            plot_approx(fn_name, r.clone(), &op1, |x| op2(Fp::from(x)).into());
            plot_err(fn_name, r.clone(), &op1, |x| op2(Fp::from(x)).into());
            if BENCH
            {
                plot_bench(fn_name, r, |x| {
                    let t0 = Instant::now();

                    for _ in 0..M
                    {
                        let _ = op1(x);
                    }

                    Instant::now().duration_since(t0).div_f64(M as f64).as_secs_f32()
                }, |x| {
                    let x = Fp::from(x);
                    let t0 = Instant::now();

                    for _ in 0..M
                    {
                        let _ = op2(x);
                    }

                    Instant::now().duration_since(t0).div_f64(M as f64).as_secs_f32()
                })
            }
        }
    }

    const M: usize = 64;
    const N: usize = 1024;
    const PLOT_TARGET: &str = "plots";
    const BENCH: bool = false;
    
    #[allow(unused)]
    pub fn plot_err<R>(
        fn_name: &str,
        range: R,
        func: impl Fn(f32) -> f32,
        approx: impl Fn(f32) -> f32
    )
    where
        R: RangeBounds<f32> + LinspaceArray<f32, N>
    {
        let x: [f32; N] = range.linspace_array();
        let y_approx = x.map(approx);

        let y = x.map(func);

        let e = y_approx.into_iter()
            .zip(y)
            .map(|(y_approx, y)| y_approx - y)
            .next_chunk()
            .unwrap();

        let plot_title: &str = &format!("{fn_name}(x) error");
        let plot_path: &str = &format!("{PLOT_TARGET}/error/{fn_name}_error.png");

        plot::plot_curves(plot_title, plot_path, [x], [e])
            .expect("Plot error");

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
    pub fn plot_approx<R>(
        fn_name: &str,
        range: R,
        func: impl Fn(f32) -> f32,
        approx: impl Fn(f32) -> f32
    )
    where
        R: RangeBounds<f32> + LinspaceArray<f32, N>
    {
        let x: [f32; N] = range.linspace_array();
        let y_approx = x.map(approx);

        let y = x.map(func);

        let plot_title: &str = &format!("{fn_name}(x)");
        let plot_path: &str = &format!("{PLOT_TARGET}/{fn_name}.png");

        plot::plot_curves(plot_title, plot_path, [x, x], [y, y_approx])
            .expect("Plot error");

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
    pub fn plot_bench<R>(
        fn_name: &str,
        range: R,
        func: impl Fn(f32) -> f32,
        approx: impl Fn(f32) -> f32
    )
    where
        R: RangeBounds<f32> + LinspaceArray<f32, N>
    {
        let x: [f32; N] = range.linspace_array();
        let y_approx = x.map(approx);

        let y = x.map(func);

        let plot_title: &str = &format!("{fn_name}(x) benchmark");
        let plot_path: &str = &format!("{PLOT_TARGET}/bench/{fn_name}_bench.png");

        plot::plot_curves(plot_title, plot_path, [x, x], [y, y_approx])
            .expect("Plot error");

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
    fn test_convert()
    {
        test_op1("ident", |x| x, |x| x, None, Some(-5.0..20.0))
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
}