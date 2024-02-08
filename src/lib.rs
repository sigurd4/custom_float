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
    use num_traits::Zero;

    use crate::{google::BFloat16, ibm::{HFPLong, HFPShort}, ieee754::{FpDouble, FpHalf, FpOctuple, FpQuadruple, FpSingle}, nvidia::TensorFloat, Fp};

    #[test]
    fn it_works() {
        let mut x = [1.0, 1.0, 0.0, 0.0].map(|x| HFPLong::from(x));
        let mut y = [Complex::zero(); 3];
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
    
    pub fn test_op2(op1: impl Fn(f32, f32) -> f32, op2: impl Fn(FpSingle, FpSingle) -> FpSingle)
    {
        for f0 in crate::tests::ttable()
        {
            for f1 in crate::tests::ttable()
            {
                let fp0 = FpSingle::from(f0);
                let fp1 = FpSingle::from(f1);

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
                    println!("{} ? {} == {} != {}", f0, f1, s, sp);
                }
            }
        }
    }
    
    pub fn test_op1(op1: impl Fn(f32) -> f32, op2: impl Fn(FpSingle) -> FpSingle)
    {
        for f0 in crate::tests::ttable()
        {
            let fp0 = FpSingle::from(f0);

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
                println!("{} ? == {} != {}", f0, s, sp);
            }
        }
    }

    #[test]
    fn test_convert()
    {
        for f0 in ttable()
        {
            let fp = FpSingle::from(f0);

            let f1: f32 = fp.into();

            if f0 != f1 && !(f0.is_nan() && f1.is_nan())
            {
                println!("{} != {}", f0, f1)
            }
        }
    }
}