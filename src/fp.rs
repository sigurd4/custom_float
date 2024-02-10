use std::cmp::Ordering;
use std::fmt::Binary;
use std::mem::size_of;
use std::num::FpCategory;

use array_math::{ArrayMath, ArrayOps};
use num_traits::{CheckedNeg, CheckedShl, CheckedShr, ConstOne, ConstZero, FloatConst, Inv, One, Signed};
use num_traits::{Float, NumCast, PrimInt, Unsigned, Zero};

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

pub trait UInt = Unsigned + PrimInt + CheckedShl + CheckedShr;
pub trait Int = Signed + PrimInt + CheckedShl + CheckedShr + CheckedNeg;

/// A custom floating point type.
/// 
/// Bit layout is as follows:
/// ```
/// No data: | Sign: | Exponent:  | Integer:  | Fractional: |
/// <  ..  > | < 1 > | <EXP_SIZE> | <INT_BIT> | <FRAC_SIZE> |
/// ```
#[derive(Clone, Copy)]
pub struct Fp<U: UInt, const EXP_SIZE: usize, const INT_BIT: bool, const FRAC_SIZE: usize>(U)
where
    [(); bitsize_of::<U>() - EXP_SIZE - INT_BIT as usize - FRAC_SIZE - 1]:;

impl<U: UInt, const EXP_SIZE: usize, const INT_BIT: bool, const FRAC_SIZE: usize> Fp<U, EXP_SIZE, INT_BIT, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - INT_BIT as usize - FRAC_SIZE - 1]:,
    [(); bitsize_of::<U>() - EXP_SIZE - false as usize - FRAC_SIZE - 1]:
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
        [(); bitsize_of::<V>() - E - false as usize - F - 1]:
    {
        if EXP_SIZE == E && INT_BIT == I
        {
            if let Some(b) = if bitsize_of::<U>() >= bitsize_of::<V>()
            {
                <U as NumCast>::from(fp.to_bits())
                    .map(|b| if FRAC_SIZE >= F
                    {
                        b << FRAC_SIZE - F
                    }
                    else
                    {
                        b >> F - FRAC_SIZE
                    })
            }
            else
            {
                let b = if FRAC_SIZE >= F
                {
                    fp.to_bits() << FRAC_SIZE - F
                }
                else
                {
                    fp.to_bits() >> F - FRAC_SIZE
                };
                <U as NumCast>::from(b)
            }
            {
                return Self::from_bits(b)
            }
        }

        let s = fp.sign_bit();

        if fp.is_nan()
        {
            if fp.is_snan()
            {
                return Self::snan()
            }
            return Self::qnan()
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
        let mut f = if df >= 0
        {
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
            Fp::from_bits(s_bit + (f >> 1usize))
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
                return if s.is_one() {Self::neg_infinity()} else {Self::infinity()}
            }

            Fp::from_bits(s_bit + f + (e << Self::EXP_POS))
        }
    }

    pub fn from_uint<I: UInt>(from: I) -> Self
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
            Fp::from_bits(f >> 1usize)
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

            Fp::from_bits(f + (e << Self::EXP_POS))
        }
    }
    
    pub fn from_int<I: Int>(from: I) -> Self
    {
        if from == I::min_value()
        {
            return Self::from_int(from + I::one()) - Self::one()
        }

        let s = from < I::zero();
        let mut e = Self::exp_bias() + NumCast::from(FRAC_SIZE).unwrap();
        let mut f = if bitsize_of::<I>() - 1 > bitsize_of::<U>()
        {
            let mut f = from.abs();

            while f >= I::one() << FRAC_SIZE + 1
            {
                e = match e.checked_add(&U::one())
                {
                    Some(e) => e,
                    None => return if s {Self::neg_infinity()} else {Self::infinity()}
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
            let mut f = <U as NumCast>::from(from.abs()).unwrap();
            
            while f >= U::one() << FRAC_SIZE + 1
            {
                e = match e.checked_add(&U::one())
                {
                    Some(e) => e,
                    None => return if s {Self::neg_infinity()} else {Self::infinity()}
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

        let s_bit = if s {U::one() << Self::SIGN_POS} else {U::zero()};
        if e.is_zero() // subnormal
        {
            Fp::from_bits(s_bit + (f >> 1usize))
        }
        else
        {
            if e > U::one() << EXP_SIZE
            {
                return if s {Self::neg_infinity()} else {Self::infinity()}
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

            Fp::from_bits(s_bit + f + (e << Self::EXP_POS))
        }
    }

    pub fn to_uint<I: UInt>(self) -> Option<I>
    {
        if self.is_zero()
        {
            return Some(I::zero())
        }
        let s = !self.sign_bit().is_zero();
        if !self.is_finite()
        {
            return None
        }
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

        let bias = Self::exp_bias();
        if e < bias
        {
            return Some(I::zero())
        }

        let shift = <u32 as NumCast>::from(e - bias).unwrap();
        
        let n = match <I as NumCast>::from(f)
        {
            Some(n) => if shift as usize > FRAC_SIZE
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
            },
            None => match <I as NumCast>::from(if shift as usize > FRAC_SIZE
            {
                match f.checked_shl(shift - FRAC_SIZE as u32)
                {
                    Some(n) => n,
                    None => return None
                }
            }
            else
            {
                f >> FRAC_SIZE - shift as usize
            })
            {
                Some(n) => n,
                None => return None
            }
        };;

        Some(n)
    }
    
    pub fn to_uint_wrapping<I: UInt>(mut self) -> I
    {
        if self.is_zero()
        {
            return I::zero()
        }
        if self.is_nan()
        {
            return I::zero()
        }
        if self.is_infinite()
        {
            return if self.is_sign_negative() {I::zero()} else {I::max_value()}
        }

        let max = Self::from_uint(I::max_value()) + Self::one();
        let min = Self::from_uint(I::min_value());
        self = max - ((max - min - ((self - min) % (max - min))) % (max - min));

        let s = !self.sign_bit().is_zero();
        if s
        {
            return I::zero()
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

        let bias = Self::exp_bias();
        if e < bias
        {
            return I::zero()
        }

        let shift = <u32 as NumCast>::from(e - bias).unwrap();
        
        let n = match <I as NumCast>::from(f)
        {
            Some(n) => if shift as usize > FRAC_SIZE
            {
                (n % (I::one() << bitsize_of::<I>() - (shift as usize - FRAC_SIZE))) << shift as usize - FRAC_SIZE
            }
            else
            {
                n >> FRAC_SIZE - shift as usize
            },
            None => <I as NumCast>::from(if shift as usize > FRAC_SIZE
            {
                (f % (U::one() << bitsize_of::<I>() - (shift as usize - FRAC_SIZE))) << shift as usize - FRAC_SIZE
            }
            else
            {
                f >> FRAC_SIZE - shift as usize
            }).unwrap()
        };

        n
    }
    
    pub fn to_int<I: Int>(self) -> Option<I>
    {
        if self.is_zero()
        {
            return Some(I::zero())
        }
        let s = !self.sign_bit().is_zero();
        if !self.is_finite()
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
        
        let bias = Self::exp_bias();
        if e < bias
        {
            return Some(I::zero())
        }

        let shift = <u32 as NumCast>::from(e - bias).unwrap();

        let n = match <I as NumCast>::from(f)
        {
            Some(n) => if shift as usize > FRAC_SIZE
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
            },
            None => match <I as NumCast>::from(if shift as usize > FRAC_SIZE
            {
                match f.checked_shl(shift - FRAC_SIZE as u32)
                {
                    Some(n) => n,
                    None => return None
                }
            }
            else
            {
                f >> FRAC_SIZE - shift as usize
            })
            {
                Some(n) => n,
                None => return None
            }
        };

        if s
        {
            n.checked_neg()
        }
        else
        {
            Some(n)   
        }
    }
    
    pub fn to_int_wrapping<I: Int>(mut self) -> I
    {
        if self.is_zero()
        {
            return I::zero()
        }
        if self.is_nan()
        {
            return I::zero()
        }
        if self.is_infinite()
        {
            return if self.is_sign_negative() {I::min_value()} else {I::max_value()}
        }

        let max = Self::from_int(I::max_value()) + Self::one();
        let min = Self::from_int(I::min_value());
        self = max - ((max - min - ((self - min) % (max - min))) % (max - min));

        let s = !self.sign_bit().is_zero();
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
        
        let bias = Self::exp_bias();
        if e < bias
        {
            return I::zero()
        }

        let shift = <u32 as NumCast>::from(e - bias).unwrap();

        let n = match <I as NumCast>::from(f)
        {
            Some(n) => if shift as usize > FRAC_SIZE
            {
                (n % (I::one() << bitsize_of::<I>() - 1 - (shift as usize - FRAC_SIZE))) << shift as usize - FRAC_SIZE
            }
            else
            {
                n >> FRAC_SIZE - shift as usize
            },
            None => <I as NumCast>::from(if shift as usize > FRAC_SIZE
            {
                (f % (U::one() << bitsize_of::<I>() - 1 - (shift as usize - FRAC_SIZE))) << shift as usize - FRAC_SIZE
            }
            else
            {
                f >> FRAC_SIZE - shift as usize
            }).unwrap()
        };

        if s
        {
            match n.checked_neg()
            {
                Some(n) => n,
                None => n
            }
        }
        else
        {
            n
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
    
    pub fn nan() -> Self
    {
        Self::snan()
    }
    
    pub fn qnan() -> Self
    {
        Self::from_bits(U::max_value() >> bitsize_of::<U>() - Self::SIGN_POS)
    }
    
    pub fn snan() -> Self
    {
        if Self::INT_POS + Self::INT_SIZE < 1
        {
            return Self::qnan()
        }
        Self::from_bits((U::max_value() >> bitsize_of::<U>() - Self::SIGN_POS) - (U::one() << Self::INT_POS + Self::INT_SIZE - 1))
    }

    pub fn is_snan(self) -> bool
    {
        if Self::INT_POS + Self::INT_SIZE < 1 || !self.is_nan()
        {
            return false
        }
        (self.to_bits() >> Self::INT_POS + Self::INT_SIZE - 1) & U::one() == U::zero()
    }

    pub fn infinity() -> Self
    {
        Self::from_bits((U::max_value() >> bitsize_of::<U>() - EXP_SIZE) << Self::EXP_POS)
    }

    pub fn neg_infinity() -> Self
    {
        -Self::infinity()
    }

    pub fn neg_zero() -> Self
    {
        -Self::zero()
    }

    pub fn min_value() -> Self
    {
        -Self::max_value()
    }

    pub fn min_positive_value() -> Self
    {
        Self::from_bits(U::one())
    }

    pub fn max_value() -> Self
    {
        Self::from_bits(((U::max_value() >> bitsize_of::<U>() - EXP_SIZE) << Self::EXP_POS) - U::one())
    }

    pub fn is_nan(self) -> bool
    {
        !self.is_finite() && !self.frac_bits().is_zero() && (!INT_BIT || !self.int_bit().is_zero())
    }

    pub fn is_infinite(self) -> bool
    {
        !self.is_finite() && self.frac_bits().is_zero() && (!INT_BIT || self.int_bit().is_zero())
    }

    pub fn is_finite(self) -> bool
    {
        self.exp_bits() != (U::max_value() >> bitsize_of::<U>() - EXP_SIZE)
    }

    pub fn is_normal(self) -> bool
    {
        if !self.is_finite()
        {
            return false
        }
        let e = self.exp_bits();
        if INT_BIT
        {
            return true
        }
        e != U::zero() && e != (U::max_value() >> bitsize_of::<U>() - EXP_SIZE)
    }

    pub fn classify(self) -> FpCategory
    {
        let e = self.exp_bits();
        if self.is_zero()
        {
            return FpCategory::Zero
        }
        if e == U::zero() && !INT_BIT
        {
            return FpCategory::Subnormal
        }
        if e == (U::max_value() >> bitsize_of::<U>() - EXP_SIZE)
        {
            if self.frac_bits().is_zero()
            {
                return FpCategory::Infinite
            }
            return FpCategory::Nan
        }
        FpCategory::Normal
    }

    pub fn floor(self) -> Self
    {
        if self.is_zero()
        {
            return Self::zero()
        }
        let s = !self.sign_bit().is_zero();
        let e = self.exp_bits();
        let mut f = self.frac_bits();
        let bias = Self::exp_bias();
        if e < bias
        {
            return if s {-Self::one()} else {Self::zero()}
        }
        let mask = match ((U::one() << FRAC_SIZE) - U::one()).checked_shr(<u32 as NumCast>::from(e - bias).unwrap())
        {
            Some(mask) => mask,
            None => U::zero()
        };
        let i = s && !(f & mask).is_zero();
        f = f & (!mask);
        let s_bit = if s {U::one() << Self::SIGN_POS} else {U::zero()};
        let n = Self::from_bits(s_bit + (e << Self::EXP_POS) + f + if INT_BIT {self.int_bit() << Self::INT_POS} else {U::zero()});
        
        if i
        {
            return n - Self::one()
        }
        n
    }

    pub fn ceil(self) -> Self
    {
        if self.is_zero()
        {
            return Self::zero()
        }
        let s = !self.sign_bit().is_zero();
        let e = self.exp_bits();
        let mut f = self.frac_bits();
        let bias = Self::exp_bias();
        if e < bias
        {
            return if s {Self::zero()} else {Self::one()}
        }
        let mask = match ((U::one() << FRAC_SIZE) - U::one()).checked_shr(<u32 as NumCast>::from(e - bias).unwrap())
        {
            Some(mask) => mask,
            None => U::zero()
        };
        let i = !s && !(f & mask).is_zero();
        f = f & (!mask);
        let s_bit = if s {U::one() << Self::SIGN_POS} else {U::zero()};
        let n = Self::from_bits(s_bit + (e << Self::EXP_POS) + f + if INT_BIT {self.int_bit() << Self::INT_POS} else {U::zero()});
        
        if i
        {
            return n + Self::one()
        }
        n
    }

    pub fn round(self) -> Self
    {
        let c = self.ceil();
        let f = self.floor();
        if if self >= Self::zero() {c - self <= self - f} else {c - self < self - f}
        {
            c
        }
        else
        {
            f
        }
    }

    pub fn trunc(self) -> Self
    {
        let s = !self.sign_bit().is_zero();
        if self.is_zero()
        {
            return if s {-Self::zero()} else {Self::zero()}
        }
        let e = self.exp_bits();
        let mut f = self.frac_bits();
        let bias = Self::exp_bias();
        if e < bias
        {
            return if s {-Self::zero()} else {Self::zero()}
        }
        let mask = match ((U::one() << FRAC_SIZE) - U::one()).checked_shr(<u32 as NumCast>::from(e - bias).unwrap())
        {
            Some(mask) => mask,
            None => U::zero()
        };
        f = f & (!mask);
        let s_bit = if s {U::one() << Self::SIGN_POS} else {U::zero()};
        return Self::from_bits(s_bit + (e << Self::EXP_POS) + f + if INT_BIT {self.int_bit() << Self::INT_POS} else {U::zero()})
    }

    pub fn fract(self) -> Self
    {
        self - self.trunc()
    }

    pub fn abs(self) -> Self
    {
        Self::from_bits(self.to_bits() & (U::max_value() >> bitsize_of::<U>() - Self::SIGN_POS))
    }

    pub fn signum(self) -> Self
    {
        if self.is_nan()
        {
            return self
        }
        let s = !self.sign_bit().is_zero();
        if s {-Self::one()} else {Self::one()}
    }

    pub fn is_sign_positive(self) -> bool
    {
        self.sign_bit().is_zero()
    }

    pub fn is_sign_negative(self) -> bool
    {
        !self.sign_bit().is_zero()
    }

    pub fn mul_add(self, a: Self, b: Self) -> Self
    {
        (self*a) + b
    }

    pub fn recip(self) -> Self
    {
        Self::one()/self
    }

    pub fn powi<I: Int>(self, mut n: I) -> Self
    {
        let mut r = Self::one();
        
        let mut x = if n < I::zero()
        {
            n = -n;
            Self::one() / self
        }
        else
        {
            self
        };
        
        if (n & I::one()) == I::one()
        {
            r = x;
        }

        loop
        {
            n = n >> 1usize;
            if n == I::zero()
            {
                break;
            }
            x = x*x;
            if (n & I::one()) == I::one()
            {
                r = r*x;
            }
        }

        r
    }
    
    pub fn powu<I: UInt>(self, mut n: I) -> Self
    {
        let mut r = Self::one();
        
        let mut x = self;
        
        if (n & I::one()) == I::one()
        {
            r = x;
        }

        loop
        {
            n = n >> 1usize;
            if n == I::zero()
            {
                break;
            }
            x = x*x;
            if (n & I::one()) == I::one()
            {
                r = r*x;
            }
        }

        r
    }

    pub fn powf(self, n: Self) -> Self
    {
        if INT_BIT
        {
            return Self::from_fp(Fp::<U, EXP_SIZE, false, FRAC_SIZE>::from_fp(self).powf(Fp::<U, EXP_SIZE, false, FRAC_SIZE>::from_fp(n)))
        }

        let xabs = self.abs();
        let nabs = n.abs();

        let edge_x = {
            let e = (U::one() << EXP_SIZE) - U::one();
            Self::from_bits(e << Self::EXP_POS)
        };
        let edge_n = {
            let e = Self::exp_bias() + U::from(FRAC_SIZE).unwrap();
            Self::from_bits(e << Self::EXP_POS)
        };
        
        if n.is_zero()
        {
            return Self::one()
        }
        if self.is_one() && !n.is_finite()
        {
            return Self::one()
        }
        if self.is_nan() || n.is_nan()
        {
            return self + n
        }
        if n == Self::one()
        {
            return self
        }
        
        // if n is an odd integer
        let n_d = (nabs + Self::one())/<Self as From<_>>::from(2.0);
        let noi = n_d.trunc() == n_d;

        if xabs.is_zero()
        {
            if noi || n.is_infinite()
            {
                if n < Self::zero()
                {
                    return Self::one()/self
                }
                return self
            }
            if n > Self::zero()
            {
                return Self::zero()
            }
            return Self::one()/xabs
        }

        if (xabs.to_bits() >= edge_x.to_bits()) || (nabs.to_bits() >= edge_n.to_bits())
        {
            if nabs.is_infinite()
            {
                if self == -Self::one()
                {
                    return Self::one()
                }
                if xabs > Self::one()
                {
                    if n > Self::zero()
                    {
                        return n
                    }
                    return Self::zero()
                }
                if n > Self::zero()
                {
                    return Self::zero()
                }
                return nabs
            }

            if self == Self::infinity()
            {
                if n < Self::zero()
                {
                    return Self::zero()
                }
                return self
            }
            if self > Self::neg_infinity()
            {
                if !noi
                {
                    if self < Self::zero() || n > Self::zero()
                    {
                        return self.sqrt()
                    }
                    return (Self::one()/self).sqrt()
                }
            
                // clamp  -0x1.0p31 < y < 0x1.0p31
                let nabs = nabs.min(Self::from_uint(U::max_value()));
                let s = !n.sign_bit().is_zero();
                let mut i = <u128 as NumCast>::from(nabs).unwrap();
                let mut r = Self::one();
                
                let mut x = if s
                {
                    Self::one() / self
                }
                else
                {
                    self
                };
                
                if (i & 1) == 1
                {
                    r = x;
                }

                loop
                {
                    i >>= 1;
                    if 0 == i
                    {
                        break;
                    }
                    x = x*x;
                    if (i & 1) == 1
                    {
                        r = r*x;
                    }
                }

                return r;
            }

            if noi
            {
                if n > Self::zero()
                {
                    return self
                }
                return -Self::zero()
            }
            
            if n > Self::zero()
            {
                return -self
            }
            return Self::zero()
        }
        
        let ni = n.trunc() == n;

        if !ni && self < Self::zero()
        {
            return Self::snan()
        }
        
        //special case for sqrts
        if <Self as From<_>>::from(0.5) == nabs
        {
            if self < Self::zero() || n > Self::zero()
            {
                return self.sqrt()
            }
            return (Self::one()/self).sqrt()
        }

        let xabs_log2 = xabs.log2();

        let n_x_abs_log2 = n*xabs_log2;
        if n_x_abs_log2 >= Self::from_uint(U::one() << (EXP_SIZE - 1))
        {
            // ???
            return Self::from_bits(((U::one() << EXP_SIZE) - U::one()) << Self::EXP_POS)*n_x_abs_log2;
        }
        if n_x_abs_log2 <= -Self::from_uint((U::one() << (EXP_SIZE - 1)) + U::from(FRAC_SIZE - 1).unwrap())
        {
            // ???
            return -Self::from_bits(U::one() << (Self::EXP_POS - EXP_SIZE/2))*n_x_abs_log2;
        }
        
        let e = n_x_abs_log2.floor();
        let f = n_x_abs_log2 - e;
        
        const C0: f64 =  1.0 + 0.278626872016317130037181614004e-10;
        const C1: f64 = 0.693147176943623740308984004029708;
        const C2: f64 = 0.240226505817268621584559118975830;
        const C3: f64 = 0.555041568519883074165425891257052e-1;
        const C4_OVER_C8: f64 = 0.961813690023115610862381719985771e-2 / C8;
        const C5_OVER_C8: f64 = 0.133318252930790403741964203236548e-2 / C8;
        const C6_OVER_C8: f64 = 0.154016177542147239746127455226575e-3 / C8;
        const C7_OVER_C8: f64 = 0.154832722143258821052933667742417e-4 / C8;
        const C8: f64 = 0.134107709538786543922336536865157e-5;

        let mut z = Self::one();
        if !f.is_zero()
        {
            let ff = f * f;
            let mut s7 = <Self as From<_>>::from(C7_OVER_C8) * f;
            let mut s3 = <Self as From<_>>::from(C3) * f;
            let mut s5 = <Self as From<_>>::from(C5_OVER_C8) * f;
            let mut s1 = <Self as From<_>>::from(C1) * f;
            let ffff = ff * ff;
            s7 += <Self as From<_>>::from(C6_OVER_C8);
            s3 += <Self as From<_>>::from(C2);
            s5 += <Self as From<_>>::from(C4_OVER_C8);
            s1 += <Self as From<_>>::from(C0);
            s7 *= ff;
            s3 *= ff;
            s5 += ffff;
            let c8ffff = ffff * <Self as From<_>>::from(C8);
            s7 += s5;
            s3 += s1;
            s7 *= c8ffff;
            z = s3 + s7;
        }

        let e = e + Self::from_uint(Self::exp_bias());
        if e > Self::from_uint((U::one() << EXP_SIZE) - U::one())
        {
            return Self::infinity()*z
        }
        if e < Self::zero()
        {
            return Self::zero()*z
        }
        if e.is_nan()
        {
            return e
        }
        let e = U::from(e).unwrap();
        let y = Self::from_bits(e << Self::EXP_POS);

        y*z
    }

    pub fn sqrt(self) -> Self
    {
        if self.is_nan()
        {
            return self
        }
        if self.is_sign_negative()
        {
            return Self::snan()
        }

        if INT_BIT
        {
            return Self::from_fp(Fp::<U, EXP_SIZE, false, FRAC_SIZE>::from_fp(self).sqrt())
        }

        let y = Self::from_bits(
            ((Self::exp_bias() + U::one()) * (U::one() << (FRAC_SIZE - 1)))
            + (self.to_bits() >> 1usize)
            - (U::one() << (FRAC_SIZE - 1))
        );

        const NEWTON: usize = 4;
        let half = <Self as From<_>>::from(0.5);
        let mut y = y;
        for _ in 0..NEWTON
        {
            y = half*(y + self/y);
        }
        y
    }

    pub fn exp(self) -> Self
    {
        (self/Self::LN_2()).exp2()
    }

    pub fn exp2(self) -> Self
    {
        if INT_BIT
        {
            return Self::from_fp(Fp::<U, EXP_SIZE, false, FRAC_SIZE>::from_fp(self).exp2())
        }
        if self.is_nan()
        {
            return self
        }
        if self.is_infinite()
        {
            if self.is_sign_negative()
            {
                return Self::zero()
            }
            return self
        }
        let neg = self.is_sign_negative();
        let x = self.abs().max(-Self::from_uint(Self::exp_bias() - U::one()));
        let w = x.floor();
        let z = x - w;
        let approx = <Self as From<_>>::from(-5.7259425)
        + <Self as From<_>>::from(27.7280233)/(<Self as From<_>>::from(4.84252568) - z)
        - <Self as From<_>>::from(1.49012907)*z;
        let y = Self::from_uint(U::one() << Self::EXP_POS)*(x + Self::from_uint(Self::exp_bias()) + approx);
        if y.is_nan()
        {
            return y
        }
        if y >= Self::from_uint(U::max_value())
        {
            if neg
            {
                return Self::zero()
            }
            return Self::infinity()
        }
        if y <= Self::one()
        {
            if neg
            {
                return Self::infinity()
            }
            return Self::zero();
        }
        let y = y.to_uint().unwrap();
        if y < U::one() << Self::EXP_POS
        {
            return Self::one()
        }
        let y = Self::from_bits(y);
        if neg
        {
            return y.recip();
        }
        y
    }

    pub fn ln(self) -> Self
    {
        self.log2()*Self::LN_2()
    }

    pub fn log(self, base: Self) -> Self
    {
        self.ln()/base.ln()
    }

    pub fn log2(self) -> Self
    {
        if INT_BIT
        {
            return Self::from_fp(Fp::<U, EXP_SIZE, false, FRAC_SIZE>::from_fp(self).log2())
        }
        if self.is_nan()
        {
            return self
        }
        if self < Self::zero()
        {
            return Self::snan()
        }
        if self.is_zero()
        {
            return Self::neg_infinity()
        }
        if self.is_infinite()
        {
            return Self::infinity()
        }
        let b = if !self.exp_bits().is_zero()
        {
            U::zero()
        }
        else
        {
            U::from(self.frac_bits().leading_zeros() as usize - (bitsize_of::<U>() - FRAC_SIZE)).unwrap()
        };
        let mut y = Self::from_uint(self.exp_bits())
            - Self::from_uint(Self::exp_bias() + U::one())
            - Self::from_uint(b);
        let mut u = self.to_bits();
        u = u & !(((U::one() << EXP_SIZE) - U::one()) << Self::EXP_POS);
        u = u + (Self::exp_bias() << Self::EXP_POS);
        let u = Fp::<U, EXP_SIZE, false, FRAC_SIZE>::from_bits(U::from(u).unwrap());
        let u = Self::from_fp(u);
        y += (<Self as From<_>>::from(-0.34484843) * u + <Self as From<_>>::from(2.02466578)) * u  - <Self as From<_>>::from(0.67487759); 
        return y;
    }

    pub fn log10(self) -> Self
    {
        self.ln()/Self::LN_10()
    }

    pub fn max(self, other: Self) -> Self
    {
        match self.partial_cmp(&other)
        {
            Some(cmp) => match cmp
            {
                Ordering::Less => other,
                _ => self,
            }
            None => self + other
        }
    }

    pub fn min(self, other: Self) -> Self
    {
        match self.partial_cmp(&other)
        {
            Some(cmp) => match cmp
            {
                Ordering::Greater => other,
                _ => self,
            }
            None => self + other
        }
    }

    pub fn abs_sub(self, other: Self) -> Self
    {
        (self - other).abs()
    }

    pub fn cbrt(self) -> Self
    {
        self.powf(Self::one()/<Self as From<_>>::from(3.0))
    }

    pub fn hypot(self, other: Self) -> Self
    {
        (self*self + other*other).sqrt()
    }

    pub fn sin(self) -> Self
    {
        const N: usize = 6;
        const C: [f64; N] = [
            1.276278962,
            -0.285261569,
            0.009118016,
            -0.000136587,
            0.000001185,
            -0.000000007
        ];
        let t: [[Self; N]; N] = ArrayOps::fill(
            |n| <[_; N]>::chebyshev_polynomial(1, n).unwrap()
        );
        let p: [Self; N] = t.zip(C)
            .map2(|(t, c)| t.map2(|tn| <Self as From<_>>::from(c)*tn))
            .reduce(|a, b| a.zip(b).map2(|(a, b)| a + b))
            .unwrap_or_default();

        let mut w = self*Self::FRAC_2_PI();
        let mut i = 0;
        while i < 4
        {
            w -= Self::one();
            if i % 2 == 0 && w < Self::zero()
            {
                w = -w;
            }
            w %= <Self as From<_>>::from(4.0);
            i += 1;
        }
        let two = <Self as From<_>>::from(2.0);
        let w = if w > Self::one() {two - w} else if w < -Self::one() {-two - w} else {w};

        let z = two*w*w - Self::one();

        p.polynomial(z)*w
    }

    pub fn cos(self) -> Self
    {
        const N: usize = 6;
        const C: [f64; N] = [
            1.276278962,
            -0.285261569,
            0.009118016,
            -0.000136587,
            0.000001185,
            -0.000000007
        ];
        let t: [[Self; N]; N] = ArrayOps::fill(
            |n| <[_; N]>::chebyshev_polynomial(1, n).unwrap()
        );
        let p: [Self; N] = t.zip(C)
            .map2(|(t, c)| t.map2(|tn| <Self as From<_>>::from(c)*tn))
            .reduce(|a, b| a.zip(b).map2(|(a, b)| a + b))
            .unwrap_or_default();

        let mut w = self*Self::FRAC_2_PI() + Self::one();
        let mut i = 0;
        while i < 4
        {
            w -= Self::one();
            if i % 2 == 0 && w < Self::zero()
            {
                w = -w;
            }
            w %= <Self as From<_>>::from(4.0);
            i += 1;
        }
        let two = <Self as From<_>>::from(2.0);
        let w = if w > Self::one() {two - w} else if w < -Self::one() {-two - w} else {w};

        let z = two*w*w - Self::one();

        p.polynomial(z)*w
    }

    pub fn tan(self) -> Self
    {
        let (sin, cos) = self.sin_cos();

        sin/cos
    }

    pub fn asin(self) -> Self
    {
        if self.is_nan()
        {
            return self
        }
        if self == Self::one()
        {
            return Self::FRAC_PI_2()
        }
        if self == -Self::one()
        {
            return -Self::FRAC_PI_2()
        }
        if self > Self::one() || self < -Self::one()
        {
            return Self::snan()
        }
        (self/(Self::one() - self*self).sqrt()).atan()
    }

    pub fn acos(self) -> Self
    {
        if self.is_nan()
        {
            return self
        }
        if self == Self::one()
        {
            return Self::zero()
        }
        if self == -Self::one()
        {
            return Self::PI()
        }
        if self > Self::one() || self < -Self::one()
        {
            return Self::snan()
        }
        ((Self::one() - self*self).sqrt()/self).atan()
    }

    pub fn atan(self) -> Self
    {
        if self.is_nan()
        {
            return self
        }
        if self.is_infinite()
        {
            return Self::FRAC_PI_2().copysign(self)
        }
        const TAYLOR: usize = 8;
        if self.abs() < Self::one()
        {
            let mut z = self;
            let mut y = z;

            for k in 1..TAYLOR
            {
                z *= -self*self;
                y += z/Self::from_uint(1 + 2*k)
            }

            y
        }
        else
        {
            let mut z = Self::one()/self;
            let mut y = Self::FRAC_PI_2() - z;

            for k in 1..TAYLOR
            {
                z /= -self*self;
                y -= z/Self::from_uint(1 + 2*k)
            }

            y
        }
    }

    pub fn atan2(self, other: Self) -> Self
    {
        if other.is_zero()
        {
            return if self.is_zero()
            {
                Self::zero()
            }
            else if self.is_sign_negative()
            {
                -Self::FRAC_PI_2()
            }
            else
            {
                Self::FRAC_PI_2()
            }
        }
        let atan = (other/self).atan();
        if other.is_sign_positive()
        {
            return atan
        }
        atan + if self.is_sign_positive()
        {
            Self::PI()
        }
        else
        {
            -Self::PI()
        }
    }

    pub fn sin_cos(self) -> (Self, Self)
    {
        const N: usize = 6;
        const C: [f64; N] = [
            1.276278962,
            -0.285261569,
            0.009118016,
            -0.000136587,
            0.000001185,
            -0.000000007
        ];
        let t: [[Self; N]; N] = ArrayOps::fill(
            |n| <[_; N]>::chebyshev_polynomial(1, n).unwrap()
        );
        let p: [Self; N] = t.zip(C)
            .map2(|(t, c)| t.map2(|tn| <Self as From<_>>::from(c)*tn))
            .reduce(|a, b| a.zip(b).map2(|(a, b)| a + b))
            .unwrap_or_default();

        let sin = {
            let mut w = self*Self::FRAC_2_PI();
            let mut i = 0;
            while i < 4
            {
                w -= Self::one();
                if i % 2 == 0 && w < Self::zero()
                {
                    w = -w;
                }
                w %= <Self as From<_>>::from(4.0);
                i += 1;
            }
            let two = <Self as From<_>>::from(2.0);
            let w = if w > Self::one() {two - w} else if w < -Self::one() {-two - w} else {w};
    
            let z = two*w*w - Self::one();
    
            p.polynomial(z)*w
        };

        let cos = {
            let mut w = self*Self::FRAC_2_PI() + Self::one();
            let mut i = 0;
            while i < 4
            {
                w -= Self::one();
                if i % 2 == 0 && w < Self::zero()
                {
                    w = -w;
                }
                w %= <Self as From<_>>::from(4.0);
                i += 1;
            }
            let two = <Self as From<_>>::from(2.0);
            let w = if w > Self::one() {two - w} else if w < -Self::one() {-two - w} else {w};

            let z = two*w*w - Self::one();

            p.polynomial(z)*w
        };

        (sin, cos)
    }

    pub fn exp_m1(self) -> Self
    {
        self.exp() - Self::one()
    }

    pub fn ln_1p(self) -> Self
    {
        (self + Self::one()).ln()
    }

    pub fn sinh(self) -> Self
    {
        if !self.is_finite()
        {
            return self
        }

        let ex = self.exp();
        let emx = (-self).exp();

        (ex - emx)*<Self as From<_>>::from(0.5)
    }

    pub fn cosh(self) -> Self
    {
        if !self.is_finite()
        {
            return self.abs()
        }

        let ex = self.exp();
        let emx = (-self).exp();

        (ex + emx)*<Self as From<_>>::from(0.5)
    }

    pub fn tanh(self) -> Self
    {
        if self.is_nan()
        {
            return self
        }

        let x2 = <Self as From<_>>::from(2.0)*self;
        let ex2 = x2.exp();
        let ex2p1 = ex2 + Self::one();
        let ex2m1 = ex2 - Self::one();

        if ex2m1.is_infinite()
        {
            return Self::one().copysign(self)
        }
        
        ex2m1/ex2p1
    }

    pub fn asinh(self) -> Self
    {
        if !self.is_finite()
        {
            return self
        }
        (self + (self*self + Self::one()).sqrt()).ln()
    }

    pub fn acosh(self) -> Self
    {
        if self.is_nan()
        {
            return self
        }
        if self < Self::one()
        {
            return Self::snan()
        }
        if self.is_infinite()
        {
            return self
        }
        (self + (self*self - Self::one()).sqrt()).ln()
    }

    pub fn atanh(self) -> Self
    {
        if self.is_nan()
        {
            return self
        }
        if self == Self::one()
        {
            return Self::infinity()
        }
        if self == -Self::one()
        {
            return Self::neg_infinity()
        }
        if self > Self::one() || self < -Self::one()
        {
            return Self::snan()
        }
        if self > Self::zero()
        {
            <Self as From<_>>::from(0.5)*((Self::one() + self)/(Self::one() - self)).ln()
        }
        else
        {
            -<Self as From<_>>::from(0.5)*((Self::one() - self)/(Self::one() + self)).ln()
        }
    }

    pub fn epsilon() -> Self
    {
        let frac_size = U::from(FRAC_SIZE).unwrap();
        let bias = Self::exp_bias();
        if bias <= frac_size
        {
            return Self::from_bits(U::one() << <usize as NumCast>::from(bias).unwrap() - 1)
        }
        Self::from_bits((bias - frac_size) << Self::EXP_POS)
    }

    pub fn copysign(self, sign: Self) -> Self
    {
        let s = sign.sign_bit();
        let e = self.exp_bits();
        let f = self.frac_bits();
        if INT_BIT
        {
            let i = self.int_bit();
            return Self::from_bits((s << Self::SIGN_POS) + (e << Self::EXP_POS) + (i << Self::INT_POS) + (f << Self::FRAC_POS))
        }
        Self::from_bits((s << Self::SIGN_POS) + (e << Self::EXP_POS) + (f << Self::FRAC_POS))
    }

    pub fn is_subnormal(self) -> bool
    {
        !INT_BIT && self.exp_bits() == U::zero() && !self.is_zero()
    }

    pub fn to_degrees(self) -> Self
    {
        self*(<Self as From<_>>::from(180.0)/Self::PI())
    }
    pub fn to_radians(self) -> Self
    {
        self*(Self::PI()/<Self as From<_>>::from(180.0))
    }
    
    pub fn zero() -> Self
    {
        Self::from_bits(U::zero())
    }

    pub fn is_zero(self) -> bool
    {
        (self.to_bits() << Self::BIT_SIZE - Self::SIGN_POS).is_zero()
    }

    fn set_zero(&mut self)
    {
        *self = Self::zero()
    }
    
    pub fn one() -> Self
    {
        let bias = Self::exp_bias();
        if INT_BIT
        {
            return Self::from_bits((bias << Self::EXP_POS) + (U::one() << Self::INT_POS))
        }
        Self::from_bits(bias << Self::EXP_POS)
    }

    pub fn is_one(self) -> bool
    {
        self == Self::one()
    }

    fn set_one(&mut self)
    {
        *self = Self::one()
    }
    
    fn div_euclid(self, rhs: Self) -> Self
    {
        let q = (self / rhs).trunc();
        if self % rhs < Self::zero()
        {
            return if rhs > Self::zero() {q - Self::one()} else {q + Self::one()};
        }
        q
    }

    fn rem_euclid(self, rhs: Self) -> Self
    {
        let r = self % rhs;
        if r < Self::zero() {r + rhs.abs()} else {r}
    }
    
    fn total_cmp(self, other: Self) -> std::cmp::Ordering
    {
        let s1 = self.sign_bit();
        let s2 = self.sign_bit();

        if s1 != s2
        {
            return s2.cmp(&s1)
        }

        let s = !s1.is_zero();
        
        let mut left = self.to_bits();
        let mut right = other.to_bits();

        left = left ^ (((left >> (bitsize_of::<U>() - 1))) >> 1usize);
        right = right ^ (((right >> (bitsize_of::<U>() - 1))) >> 1usize);

        if s {right.cmp(&left)} else {left.cmp(&right)}
    }
}

impl<U: UInt, const EXP_SIZE: usize, const INT_BIT: bool, const FRAC_SIZE: usize> Fp<U, EXP_SIZE, INT_BIT, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - INT_BIT as usize - FRAC_SIZE - 1]:,
    [(); bitsize_of::<U>() - EXP_SIZE - false as usize - FRAC_SIZE - 1]:,
    U: ConstZero
{
    pub const ZERO: Self = Self::from_bits(U::ZERO);
}

impl<U: UInt, const EXP_SIZE: usize, const INT_BIT: bool, const FRAC_SIZE: usize> Fp<U, EXP_SIZE, INT_BIT, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - INT_BIT as usize - FRAC_SIZE - 1]:,
    [(); bitsize_of::<U>() - EXP_SIZE - false as usize - FRAC_SIZE - 1]:,
    U: ConstOne
{
    pub const MIN_POSITIVE: Self = Self::from_bits(U::ONE);
}

pub const fn bitsize_of<T>() -> usize
{
    size_of::<T>()*8
}