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
        util
    }
);

pub trait UInt = Unsigned + PrimInt + CheckedShl + CheckedShr + Debug;
pub trait Int = Signed + PrimInt + CheckedShl + CheckedShr + CheckedNeg;

#[cfg(test)]
mod tests
{
    #![allow(unused)]

    use array_math::ArrayMath;
    use num::Complex;
    use num_traits::{Float, One, ToPrimitive, Zero};

    use crate::{amd::Fp24, g_711::FpG711, google::{Bf16, Bf32, Bf8}, ibm::{HFpLong, HFpShort}, ieee754::{DecDouble, FpDouble, FpHalf, FpQuadruple, FpSingle}, intel::Fp80, nvidia::Tf19, Fp};

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
    
    pub fn test_op2(op1: impl Fn(f32, f32) -> f32, op2: impl Fn(F, F) -> F, d: Option<f32>)
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
    
    pub fn test_op1(op1: impl Fn(f32) -> f32, op2: impl Fn(F) -> F, d: Option<f32>)
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