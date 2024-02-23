#![allow(incomplete_features)]
#![feature(const_trait_impl)]
#![feature(generic_const_exprs)]
#![feature(trait_alias)]
#![feature(const_fn_floating_point_arithmetic)]
#![feature(portable_simd)]
#![feature(associated_type_bounds)]
#![feature(round_ties_even)]
#![feature(float_gamma)]
#![feature(lazy_cell)]

use std::fmt::Debug;

use num_traits::{CheckedNeg, CheckedShl, CheckedShr, PrimInt, Signed, Unsigned};

moddef::moddef!(
    pub mod {
        amd,
        g_711,
        google,
        ibm,
        ieee754,
        intel,
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

    use std::ops::{RangeBounds, Range};

    use array_math::{ArrayMath, ArrayOps};
    use linspace::LinspaceArray;
    use num::Complex;
    use num_traits::{Float, One, ToPrimitive, Zero};

    use crate::{amd::Fp24, g_711::FpG711, google::{Bf16, Bf32, Bf8}, ibm::{HFpLong, HFpShort}, ieee754::{DecDouble, FpDouble, FpHalf, FpQuadruple, FpSingle}, intel::Fp80, nvidia::Tf19, plot, Fp};

    pub type F = FpSingle;

    #[test]
    fn it_works()
    {
        let mut x = [0.0; 8].map(|x| F::from(x));
        x[0] = F::one() + F::one() + F::one();
        x[1] = F::one() + F::one();
        x[2] = F::one();

        println!("{:?}", x);

        let mut y = [Complex::zero(); 5];
        x.real_fft(&mut y);
        x.real_ifft(&y);

        println!("{:?}", x);
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
            plot_err(fn_name, r, op1, |x| op2(Fp::from(x)).into())
        }
    }

    const N: usize = 1024;
    const PLOT_TARGET: &str = "plots";
    
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

        let e = y_approx.sub_each(y);

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
        for f0 in ttable::<f32>()
        {
            let fp = F::from(f0);

            let f1: f32 = fp.into();

            if f0 != f1 && !(f0.is_nan() && f1.is_nan())
            {
                println!("{} != {}", f0, f1)
            }
        }
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