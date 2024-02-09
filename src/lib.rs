#![allow(incomplete_features)]
#![feature(const_trait_impl)]
#![feature(generic_const_exprs)]
#![feature(trait_alias)]
#![feature(const_fn_floating_point_arithmetic)]
#![feature(portable_simd)]

moddef::moddef!(
    flat(pub) mod {
        fp
    }
);

#[cfg(test)]
mod tests
{
    use array_math::ArrayMath;
    use num::Complex;
    use num_traits::{One, Zero};

    use crate::{amd::Fp24, g_711::FpG711, google::{Bf16, Bf32, Bf8}, ibm::{HFpLong, HFpShort}, ieee754::{FpDouble, FpHalf, FpOctuple, FpQuadruple, FpSingle}, intel::Fp80, nvidia::{TensorFloat19, TensorFloat32}, Fp};

    type F = Bf32;

    #[test]
    fn it_works()
    {
        let mut x = [0.0; 8].map(|x| F::from(x));
        x[0] = One::one();
        x[1] = One::one();
        x[2] = One::one();

        println!("{:?}", x);

        let mut y = [Complex::zero(); 5];
        x.real_fft(&mut y);
        x.real_ifft(&y);

        println!("{:?}", x);
    }

    pub fn ttable() -> Vec<f32>
    {
        vec![
            3.333333,
            10.0,
            16.0,
            -2.2,
            2.2,
            1.0,
            1.0 + f32::EPSILON,
            0.0,
            0.5,
            -0.5,
            f32::NAN,
            f32::INFINITY,
            f32::NEG_INFINITY,
            f32::MIN_POSITIVE,
            f32::MIN_POSITIVE/2.0,
            f32::from_bits(1),
        ]
    }
    
    pub fn test_op2(op1: impl Fn(f32, f32) -> f32, op2: impl Fn(F, F) -> F)
    {
        for f0 in crate::tests::ttable()
        {
            for f1 in crate::tests::ttable()
            {
                let fp0 = F::from(f0);
                let fp1 = F::from(f1);

                let s = op1(f0, f1);
                let sp = op2(fp0, fp1).into();

                if s != sp && !(s.is_nan() && sp.is_nan())
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
    
    pub fn test_op1(op1: impl Fn(f32) -> f32, op2: impl Fn(F) -> F)
    {
        for f0 in crate::tests::ttable()
        {
            let fp0 = F::from(f0);

            let s = op1(f0);
            let sp = op2(fp0).into();

            if s != sp && !(s.is_nan() && sp.is_nan())
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
    fn to_int()
    {
        let x = FpSingle::from_uint(3usize);
        println!("{}", x.to_uint::<u64>().unwrap());
        let x = f32::from_bits(x.to_bits());
        println!("{}", x);
    }

    #[test]
    fn int_bit()
    {
        let x = Fp80::from(2.0);
        println!("{}", x.int_bit());
        println!("{}", x);
        println!("{:b}", x)
    }

    #[test]
    fn test_convert()
    {
        for f0 in ttable()
        {
            let fp = F::from(f0);

            let f1: f32 = fp.into();

            if f0 != f1 && !(f0.is_nan() && f1.is_nan())
            {
                println!("{} != {}", f0, f1)
            }
        }
    }
}