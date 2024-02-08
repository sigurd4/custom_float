use std::fmt::{Binary, Display};
use std::mem::size_of;

use num_traits::{CheckedShl, CheckedShr, NumAssign, ToPrimitive};
use num_traits::{ops::overflowing::OverflowingMul, Float, Num, NumCast, PrimInt, Signed, Unsigned, Zero};

moddef::moddef!(
    pub mod {
        amd,
        g_711,
        google,
        ibm,
        ieee754,
        nvidia,
        pixar
    },
    flat(pub) mod {
        binary,
        debug,
        default,
        display,
        float,
        float_const,
        float_core,
        num_,
        num_cast,
        one,
        ops,
        partial_eq,
        partial_ord,
        to_primitive,
        zero,
    }
);

pub trait UInt = Unsigned + PrimInt + OverflowingMul + Binary + CheckedShl + CheckedShr;

/// A custom floating point type.
/// 
/// Bit layout is as follows:
/// ```
/// No data: | Sign: | Exponent:  | Fractional: |
/// <  ..  > | < 1 > | <EXP_SIZE> | <FRAC_SIZE> |
/// ```
#[derive(Clone, Copy)]
pub struct Fp<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize>(U)
where
    [(); bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 1]:;

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> Fp<U, EXP_SIZE, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 1]:
{
    const BIT_SIZE: usize = EXP_SIZE + FRAC_SIZE + 1;
    const SIGN_SIZE: usize = 1;
    const SIGN_POS: usize = EXP_SIZE + FRAC_SIZE;
    const EXP_POS: usize = FRAC_SIZE;
    const FRAC_POS: usize = 0;

    pub fn from_fp<V: UInt, const E: usize, const F: usize>(fp: Fp<V, E, F>) -> Self
    where
        [(); bitsize_of::<V>() - E - F - 1]:
    {
        let s = fp.sign_bit();

        if fp.is_nan()
        {
            return Self::nan()
        }
        if fp.is_infinite()
        {
            return if s == V::one() {Self::neg_infinity()} else {Self::infinity()}
        }
        if fp.is_zero()
        {
            return if s == V::one() {Self::neg_zero()} else {Self::zero()}
        }

        let e = fp.exp_bits();
        let f = fp.frac_bits();

        let df = FRAC_SIZE as isize - F as isize;

        let bias1 = Fp::<V, E, F>::exp_bias();
        let bias2 = Self::exp_bias();
        let s = U::from(s).unwrap();
        let mut f = if df >= 0 {
            U::from(f).unwrap() << df as usize
        }
        else
        {
            NumCast::from(f >> (-df) as usize).unwrap()
        };
        if !e.is_zero() //normal
        {
            f = f + (U::one() << FRAC_SIZE);
        }
        else
        {
            f = f << 1usize;
        }
        let mut e = if Self::BIT_SIZE >= Fp::<V, E, F>::BIT_SIZE
        {
            let bias1 = U::from(bias1).unwrap();
            if bias2 >= bias1
            {
                U::from(e).unwrap() + (bias2 - bias1)
            }
            else
            {
                let e = U::from(e).unwrap();
                match e.checked_sub(&(bias1 - bias2))
                {
                    Some(e) => e,
                    None => {
                        let shift = NumCast::from(bias1 - bias2 - e - U::one()).unwrap();
                        f = match f.checked_shr(shift)
                        {
                            Some(f) => f,
                            None => U::zero()
                        };
                        U::zero()
                    }
                }
            }
        }
        else
        {
            let bias2 = V::from(bias2).unwrap();
            
            if bias2 >= bias1
            {
                U::from(e + (bias2 - bias1)).unwrap()
            }
            else
            {
                match e.checked_sub(&(bias1 - bias2))
                {
                    Some(e) => U::from(e).unwrap(),
                    None => {
                        let shift = NumCast::from(bias1 - bias2 - e - V::one()).unwrap();
                        f = match f.checked_shr(shift)
                        {
                            Some(f) => f,
                            None => U::zero()
                        };
                        U::zero()
                    }
                }
            }
        };
        
        let s_bit = s << Self::SIGN_POS;
        
        while f >= U::one() << FRAC_SIZE + 1
        {
            e = e + U::one();
            f = f >> 1usize;
        }
        while e > U::zero() && f <= U::one() << FRAC_SIZE
        {
            e = e - U::one();
            f = f << 1usize;
        }

        if e.is_zero() // subnormal
        {
            return Fp::from_bits(s_bit + (f >> 1usize))
        }
        else
        {
            if e > U::one() << EXP_SIZE
            {
                return if s.is_one() {Self::negative_infinity()} else {Self::infinity()}
            }

            f = f - (U::one() << FRAC_SIZE);

            return Fp::from_bits(s_bit + f + (e << Self::EXP_POS))
        }
    }

    pub const fn from_bits(bits: U) -> Self
    {
        Self(bits)
    }
    pub const fn to_bits(self) -> U
    {
        self.0
    }

    fn sign_bit(self) -> U
    {
        (self.to_bits() & U::max_value() >> bitsize_of::<U>() - Self::SIGN_POS - Self::SIGN_SIZE) >> Self::SIGN_POS
    }
    fn exp_bits(self) -> U
    {
        assert!(EXP_SIZE + 1 < Self::BIT_SIZE);
        (self.to_bits() & U::max_value() >> bitsize_of::<U>() - Self::EXP_POS - EXP_SIZE) >> Self::EXP_POS
    }
    fn frac_bits(self) -> U
    {
        assert!(EXP_SIZE + 1 < Self::BIT_SIZE);
        (self.to_bits() & U::max_value() >> bitsize_of::<U>() - Self::FRAC_POS - FRAC_SIZE) >> Self::FRAC_POS
    }

    fn exp_bias() -> U
    {
        U::max_value() >> bitsize_of::<U>() + 1 - EXP_SIZE
    }

    fn infinity() -> Self
    {
        Self::from_bits((U::max_value() >> bitsize_of::<U>() - EXP_SIZE) << Self::EXP_POS)
    }
    fn negative_infinity() -> Self
    {
        -Self::infinity()
    }
    fn nan() -> Self
    {
        Self::from_bits(U::max_value() >> bitsize_of::<U>() - Self::SIGN_POS)
    }

    /*#[deprecated]
    pub fn to_float<F: Float>(self) -> F
    {
        let sign = !self.sign_bit().is_zero();

        if self.is_zero()
        {
            if sign {F::neg_zero()} else {F::zero()}
        }
        else if self.is_finite()
        {
            let is_normal = self.exp_bits().is_zero();
            let frac = {
                if is_normal {F::zero()} else {F::one()} //denorm
            } + F::from(self.frac_bits()).unwrap()*F::exp2(-F::from(FRAC_SIZE).unwrap());
            let bias = U::max_value() >> bitsize_of::<U>() + 1 - EXP_SIZE;
            let exp = self.exp_bits();

            return {
                if sign {-frac} else {frac}
            }*F::exp2(if is_normal {F::one()} else {F::zero()} + {
                if exp >= bias
                {
                    F::from(exp - bias).unwrap()
                }
                else
                {
                    -F::from(bias - exp).unwrap()
                }
            })
        }
        else if self.frac_bits().is_zero()
        {
            if sign {F::neg_infinity()} else {F::infinity()}
        }
        else
        {
            F::nan()
        }
    }*/
}

pub const fn bitsize_of<T>() -> usize
where
    T: UInt
{
    size_of::<T>()*8
}