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
        intel,
        nvidia,
        pixar
    },
    mod {
        impls
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
pub struct Fp<U: UInt, const EXP_SIZE: usize, const INT_BIT: bool, const FRAC_SIZE: usize>(U)
where
    [(); bitsize_of::<U>() - EXP_SIZE - INT_BIT as usize - FRAC_SIZE - 1]:;

impl<U: UInt, const EXP_SIZE: usize, const INT_BIT: bool, const FRAC_SIZE: usize> Fp<U, EXP_SIZE, INT_BIT, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - INT_BIT as usize - FRAC_SIZE - 1]:
{
    const BIT_SIZE: usize = EXP_SIZE + Self::INT_SIZE + FRAC_SIZE + 1;
    const SIGN_SIZE: usize = 1;
    const SIGN_POS: usize = EXP_SIZE + Self::INT_SIZE + FRAC_SIZE;
    const EXP_POS: usize = Self::INT_SIZE + FRAC_SIZE;
    const INT_SIZE: usize = INT_BIT as usize;
    const INT_POS: usize = FRAC_SIZE;
    const FRAC_POS: usize = 0;

    pub fn from_fp<V: UInt, const E: usize, const I: bool, const F: usize>(fp: Fp<V, E, I, F>) -> Self
    where
        [(); bitsize_of::<V>() - E - I as usize - F - 1]:,
        [(); bitsize_of::<V>() - E - false as usize - F - 1]:,
        [(); bitsize_of::<U>() - EXP_SIZE - false as usize - FRAC_SIZE - 1]:
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

        let mut e = fp.exp_bits();
        let mut f = fp.frac_bits();

        let df = FRAC_SIZE as isize - F as isize;

        let bias1 = Fp::<V, E, I, F>::exp_bias();
        let bias2 = Self::exp_bias();
        let s = U::from(s).unwrap();
        if !e.is_zero() || INT_BIT //normal
        {
            f = f + (fp.int_bit() << Fp::<V, E, I, F>::INT_POS);
        }
        else
        {
            f = f << 1usize;
        }
        let mut f = if df >= 0 {
            U::from(f).unwrap() << df as usize
        }
        else
        {
            NumCast::from(f >> (-df) as usize).unwrap()
        };
        let mut e = if Self::BIT_SIZE >= Fp::<V, E, I, F>::BIT_SIZE
        {
            let bias1 = U::from(bias1).unwrap();
            if bias2 >= bias1
            {
                U::from(e).unwrap() + (bias2 - bias1)
            }
            else
            {
                let mut e = U::from(e).unwrap();
                loop
                {
                    match e.checked_sub(&(bias1 - bias2))
                    {
                        Some(e) => break e,
                        None => {
                            e = e + U::one();
                            f = f >> 1usize;
                        }
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
                loop
                {
                    match e.checked_sub(&(bias1 - bias2))
                    {
                        Some(e) => break U::from(e).unwrap(),
                        None => {
                            e = e + V::one();
                            f = f >> 1usize;
                        }
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
            if !INT_BIT
            {
                f = f - (U::one() << Self::INT_POS);
            }
            else
            {
                f = f >> 1usize;
                e = e + U::one();
            }
            
            if e > U::one() << EXP_SIZE
            {
                return if s.is_one() {Self::negative_infinity()} else {Self::infinity()}
            }

            return Fp::from_bits(s_bit + f + (e << Self::EXP_POS))
        }
    }

    pub fn from_uint<I: UInt>(from: I) -> Self
    where
        [(); bitsize_of::<U>() - EXP_SIZE - false as usize - FRAC_SIZE - 1]:
    {
        let mut e = Self::exp_bias() + NumCast::from(FRAC_SIZE).unwrap();
        let mut f = if bitsize_of::<I>() > bitsize_of::<U>()
        {
            let mut f = from;

            while f >= I::one() << FRAC_SIZE + 1
            {
                e = match e.checked_add(&U::one())
                {
                    Some(e) => e,
                    None => return Self::infinity()
                };
                f = f >> 1usize;
            }
            while e > U::zero() && f <= I::one() << FRAC_SIZE
            {
                e = e - U::one();
                f = f << 1usize;
            }

            NumCast::from(f).unwrap()
        }
        else
        {
            let mut f = <U as NumCast>::from(from).unwrap();
            
            while f >= U::one() << FRAC_SIZE + 1
            {
                e = match e.checked_add(&U::one())
                {
                    Some(e) => e,
                    None => return Self::infinity()
                };
                f = f >> 1usize;
            }
            while e > U::zero() && f <= U::one() << FRAC_SIZE
            {
                e = e - U::one();
                f = f << 1usize;
            }

            f
        };

        if e.is_zero() // subnormal
        {
            return Fp::from_bits(f >> 1usize)
        }
        else
        {
            if e > U::one() << EXP_SIZE
            {
                return Self::infinity()
            }

            if !INT_BIT
            {
                f = f - (U::one() << Self::INT_POS);
            }
            else
            {
                f = f >> 1usize;
                e = e + U::one();
            }

            return Fp::from_bits(f + (e << Self::EXP_POS))
        }
    }

    pub fn to_uint<I: UInt>(self) -> Option<I>
    where
        [(); bitsize_of::<U>() - EXP_SIZE - false as usize - FRAC_SIZE - 1]:
    {
        if self.is_zero()
        {
            return Some(I::zero())
        }
        let s = !self.sign_bit().is_zero();
        if s
        {
            return None
        }
        let e = self.exp_bits();
        let mut f = self.frac_bits();
        
        if !e.is_zero() || INT_BIT //normal
        {
            f = f + (self.int_bit() << Self::INT_POS);
        }
        else
        {
            f = f << 1usize;
        }

        let mut n = match <I as NumCast>::from(f)
        {
            Some(n) => n,
            None => return None
        };

        let bias = Self::exp_bias();
        if e < bias
        {
            return Some(I::zero())
        }

        let shift = <u32 as NumCast>::from(e - bias).unwrap();
        n = if shift as usize > FRAC_SIZE
        {
            match n.checked_shl(shift - FRAC_SIZE as u32)
            {
                Some(n) => n,
                None => return None
            }
        }
        else
        {
            n >> FRAC_SIZE - shift as usize
        };

        Some(n)
    }

    pub const fn from_bits(bits: U) -> Self
    {
        Self(bits)
    }
    pub const fn to_bits(self) -> U
    {
        self.0
    }

    pub fn sign_bit(self) -> U
    {
        (self.to_bits() & U::max_value() >> bitsize_of::<U>() - Self::SIGN_POS - Self::SIGN_SIZE) >> Self::SIGN_POS
    }
    pub fn exp_bits(self) -> U
    {
        (self.to_bits() & U::max_value() >> bitsize_of::<U>() - Self::EXP_POS - EXP_SIZE) >> Self::EXP_POS
    }
    pub fn int_bit(self) -> U
    {
        if !INT_BIT
        {
            return U::one()
        }
        (self.to_bits() & U::max_value() >> bitsize_of::<U>() - Self::INT_POS - Self::INT_SIZE) >> Self::INT_POS
    }
    pub fn frac_bits(self) -> U
    {
        (self.to_bits() & U::max_value() >> bitsize_of::<U>() - Self::FRAC_POS - FRAC_SIZE) >> Self::FRAC_POS
    }

    pub fn exp_bias() -> U
    {
        U::max_value() >> bitsize_of::<U>() + 1 - EXP_SIZE
    }

    pub fn infinity() -> Self
    {
        Self::from_bits((U::max_value() >> bitsize_of::<U>() - EXP_SIZE) << Self::EXP_POS)
    }
    pub fn negative_infinity() -> Self
    {
        -Self::infinity()
    }
    pub fn nan() -> Self
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