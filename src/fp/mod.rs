use std::cmp::Ordering;
use std::num::FpCategory;

use array_math::{ArrayMath, ArrayOps};
use num_traits::{ConstOne, ConstZero, FloatConst};
use num_traits::NumCast;

use crate::{util, Int, UInt};

moddef::moddef!(
    mod {
        cast,
        float,
        fmt,
        identities,
        ops,

        default,
        num,
        partial_eq,
        partial_ord,
        signed,
    }
);

/// A custom floating point type.
/// 
/// Bit layout is as follows:
/// ```txt
/// No data: | Sign: | Exponent:  | Integer:   | Fractional: |
/// <  ..  > | < 1 > | <EXP_SIZE> | <INT_SIZE> | <FRAC_SIZE> |
/// ```
#[derive(Clone, Copy)]
pub struct Fp<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize>(U)
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); EXP_BASE - 2]:;

impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:,
    [(); EXP_BASE - 2]:
{
    pub const BIT_SIZE: usize = EXP_SIZE + INT_SIZE + FRAC_SIZE + 1;
    pub const SIGN_SIZE: usize = 1;

    pub const SIGN_POS: usize = EXP_SIZE + INT_SIZE + FRAC_SIZE;
    pub const EXP_POS: usize = INT_SIZE + FRAC_SIZE;
    pub const INT_POS: usize = FRAC_SIZE;
    pub const FRAC_POS: usize = 0;

    pub const IS_INT_IMPLICIT: bool = INT_SIZE == 0;

    const MANTISSA_OP_SIZE: usize = FRAC_SIZE + INT_SIZE + Self::IS_INT_IMPLICIT as usize;
    const BASE_SIZE: usize = util::bitsize_of::<usize>() - EXP_BASE.leading_zeros() as usize - 1;

    pub fn from_fp<V: UInt, const E: usize, const I: usize, const F: usize, const B: usize>(fp: Fp<V, E, I, F, B>) -> Self
    where
        [(); util::bitsize_of::<V>() - E - I - F - 1]:,
        [(); util::bitsize_of::<V>() - E - 0 - F - 1]:,
        [(); B - 2]:
    {
        /*if EXP_SIZE == E && INT_SIZE == I && EXP_BASE == B
        {
            if let Some(b) = if util::bitsize_of::<U>() >= util::bitsize_of::<V>()
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
        }*/

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

        let mut e1 = fp.exp_bits();
        let mut f = fp.frac_bits();

        let s = U::from(s).unwrap();
        if !e1.is_zero() || !Fp::<V, E, I, F, B>::IS_INT_IMPLICIT //normal
        {
            f = f + (fp.int_bits() << Fp::<V, E, I, F, B>::INT_POS);
        }
        else
        {
            f = f << 1usize;
        }

        let df = FRAC_SIZE as isize - F as isize;

        let base1 = V::from(B).unwrap();

        let mut f = loop
        {
            match if df >= 0
            {
                U::from(f).and_then(|f| if f.leading_zeros() as usize >= df as usize
                {
                    f.checked_shl(df as u32)
                }
                else
                {
                    None
                })
            }
            else
            {
                U::from(util::rounding_div_pow(f, V::from(2).unwrap(), (-df) as usize))
            }
            {
                Some(f) => break f,
                None => {
                    e1 = e1 + V::one();
                    f = util::rounding_div(f, base1)
                }
            }
        };
        
        let bias1 = Fp::<V, E, I, F, B>::exp_bias();
        let bias2 = Self::exp_bias();

        let base1 = U::from(B).unwrap();
        let base2 = U::from(EXP_BASE).unwrap();

        let mut e = bias2;
        while e1 > bias1
        {
            e1 = e1 - V::one();
            while f.leading_zeros() as usize <= Fp::<V, E, I, F, B>::BASE_SIZE
            {
                e = e + U::one();
                f = util::rounding_div(f, base2);
            }
            f = f*base1
        }
        while e1 < bias1
        {
            e1 = e1 + V::one();
            while e > U::zero() && f.leading_zeros() as usize > Self::BASE_SIZE
            {
                e = e - U::one();
                f = f*base2;
            }
            f = util::rounding_div(f, base1)
        }
        
        let s_bit = s << Self::SIGN_POS;

        while e > U::zero() && f < U::one() << Self::MANTISSA_OP_SIZE - Self::BASE_SIZE
        {
            e = e - U::one();
            f = f*base2;
        }
        while e < U::one() << EXP_SIZE && f >= U::one() << Self::MANTISSA_OP_SIZE
        {
            e = e + U::one();
            f = util::rounding_div(f, base2);
        }

        if e.is_zero() && Self::IS_INT_IMPLICIT // subnormal
        {
            Self::from_bits(s_bit + util::rounding_div_2(f))
        }
        else
        {            
            if e >= (U::one() << EXP_SIZE) - U::one()
            {
                return if !s.is_zero() {Self::neg_infinity()} else {Self::infinity()}
            }
            
            if Self::IS_INT_IMPLICIT
            {
                f = f - (U::one() << FRAC_SIZE);
                assert!(f < (U::one() << FRAC_SIZE))
            }

            Self::from_bits(s_bit + f + (e << Self::EXP_POS))
        }
    }

    pub fn from_uint<I: UInt>(from: I) -> Self
    {
        let mut e = Self::exp_bias();
        
        let mut f = if util::bitsize_of::<I>() > util::bitsize_of::<U>()
        {
            let base = I::from(EXP_BASE).unwrap();
            let mut f = from;
                
            for _ in 0..FRAC_SIZE
            {
                f = match f.checked_shl(1)
                {
                    Some(x) => x,
                    None => {
                        let x = util::rounding_div(f, base);
                        e = match e.checked_add(&U::one())
                        {
                            Some(e) => e,
                            None => return Self::infinity()
                        };
                        x << 1usize
                    }
                };
            }

            while e > U::zero() && f < I::one() << Self::MANTISSA_OP_SIZE - Self::BASE_SIZE
            {
                e = e - U::one();
                f = f*base;
            }
            while f >= I::one() << Self::MANTISSA_OP_SIZE
            {
                e = match e.checked_add(&U::one())
                {
                    Some(e) => e,
                    None => return Self::infinity()
                };
                f = util::rounding_div(f, base);
            }

            NumCast::from(f).unwrap()
        }
        else
        {
            let base = U::from(EXP_BASE).unwrap();
            let mut f = <U as NumCast>::from(from).unwrap();
                
            for _ in 0..FRAC_SIZE
            {
                f = match f.checked_shl(1)
                {
                    Some(x) => x,
                    None => {
                        let x = util::rounding_div(f, base);
                        e = match e.checked_add(&U::one())
                        {
                            Some(e) => e,
                            None => return Self::infinity()
                        };
                        x << 1usize
                    }
                };
            }
            
            while e > U::zero() && f < U::one() << Self::MANTISSA_OP_SIZE - Self::BASE_SIZE
            {
                e = e - U::one();
                f = f*base;
            }
            while f >= U::one() << Self::MANTISSA_OP_SIZE
            {
                e = match e.checked_add(&U::one())
                {
                    Some(e) => e,
                    None => return Self::infinity()
                };
                f = util::rounding_div(f, base);
            }

            f
        };

        if e.is_zero() && Self::IS_INT_IMPLICIT // subnormal
        {
            Fp::from_bits(util::rounding_div_2(f))
        }
        else
        {
            if e >= (U::one() << EXP_SIZE) - U::one()
            {
                return Self::infinity()
            }

            if Self::IS_INT_IMPLICIT
            {
                f = f - (U::one() << Self::INT_POS);
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
        let mut e = Self::exp_bias();
        let mut f = if util::bitsize_of::<I>() - 1 > util::bitsize_of::<U>()
        {
            let base = I::from(EXP_BASE).unwrap();
            let mut f = from.abs();
            
            for _ in 0..FRAC_SIZE
            {
                f = match f.checked_shl(1)
                {
                    Some(x) => x,
                    None => {
                        let x = util::rounding_div(f, base);
                        e = match e.checked_add(&U::one())
                        {
                            Some(e) => e,
                            None => return Self::infinity()
                        };
                        x << 1usize
                    }
                };
            }

            while e > U::zero() && f < I::one() << Self::MANTISSA_OP_SIZE - Self::BASE_SIZE
            {
                e = e - U::one();
                f = f*base;
            }
            while f >= I::one() << Self::MANTISSA_OP_SIZE
            {
                e = match e.checked_add(&U::one())
                {
                    Some(e) => e,
                    None => return if s {Self::neg_infinity()} else {Self::infinity()}
                };
                f = util::rounding_div(f, base);
            }

            NumCast::from(f).unwrap()
        }
        else
        {
            let base = U::from(EXP_BASE).unwrap();
            let mut f = <U as NumCast>::from(from.abs()).unwrap();
            
            for _ in 0..FRAC_SIZE
            {
                f = match f.checked_shl(1)
                {
                    Some(x) => x,
                    None => {
                        let x = util::rounding_div(f, base);
                        e = match e.checked_add(&U::one())
                        {
                            Some(e) => e,
                            None => return Self::infinity()
                        };
                        x << 1usize
                    }
                };
            }

            while e > U::zero() && f < U::one() << Self::MANTISSA_OP_SIZE - Self::BASE_SIZE
            {
                e = e - U::one();
                f = f*base;
            }
            while f >= U::one() << Self::MANTISSA_OP_SIZE
            {
                e = match e.checked_add(&U::one())
                {
                    Some(e) => e,
                    None => return if s {Self::neg_infinity()} else {Self::infinity()}
                };
                f = util::rounding_div(f, base);
            }

            f
        };

        let s_bit = if s {U::one() << Self::SIGN_POS} else {U::zero()};
        if e.is_zero() && Self::IS_INT_IMPLICIT // subnormal
        {
            Fp::from_bits(s_bit + util::rounding_div_2(f))
        }
        else
        {
            if e >= (U::one() << EXP_SIZE) - U::one()
            {
                return if s {Self::neg_infinity()} else {Self::infinity()}
            }

            if Self::IS_INT_IMPLICIT
            {
                f = f - (U::one() << Self::INT_POS);
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
        let mut e = self.exp_bits();
        let mut f = self.frac_bits();
        
        if !e.is_zero() || !Self::IS_INT_IMPLICIT //normal
        {
            f = f + (self.int_bits() << Self::INT_POS);
        }
        else
        {
            f = f << 1usize;
        }

        let bias = Self::exp_bias();

        let n = if util::bitsize_of::<I>() > util::bitsize_of::<U>()
        {
            let mut f = I::from(f).unwrap();
            let base = I::from(EXP_BASE).unwrap();

            for _ in 0..FRAC_SIZE
            {
                while e > bias && f.leading_zeros() > 3
                {
                    e = e - U::one();
                    f = f*base;
                }
                f = util::rounding_div_2(f);
            }

            while e > bias
            {
                e = e - U::one();
                f = f*base;
            }
            while e < bias
            {
                e = e + U::one();
                f = util::rounding_div(f, base);
            }

            f
        }
        else
        {
            let base = U::from(EXP_BASE).unwrap();
            let mut f = U::from(f).unwrap();

            for _ in 0..FRAC_SIZE
            {
                while e > bias && f.leading_zeros() > 3
                {
                    e = e - U::one();
                    f = f*base;
                }
                f = util::rounding_div_2(f);
            }

            while e > bias
            {
                e = e - U::one();
                f = f*base;
            }
            while e < bias
            {
                e = e + U::one();
                f = util::rounding_div(f, base);
            }

            match I::from(f)
            {
                Some(f) => f,
                None => return None
            }
        };

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
        self = max - ((max - min - ((self - min + Self::one()) % (max - min))) % (max - min)) - Self::one();

        let s = !self.sign_bit().is_zero();
        if s
        {
            return I::zero()
        }
        let mut e = self.exp_bits();
        let mut f = self.frac_bits();
        
        if !e.is_zero() || INT_SIZE != 0 //normal
        {
            f = f + (self.int_bits() << Self::INT_POS);
        }
        else
        {
            f = f << 1usize;
        }

        let bias = Self::exp_bias();
        
        if util::bitsize_of::<I>() > util::bitsize_of::<U>()
        {
            let mut f = I::from(f).unwrap();
            let base = I::from(EXP_BASE).unwrap();

            for _ in 0..FRAC_SIZE
            {
                while e > bias && f.leading_zeros() > 3
                {
                    e = e - U::one();
                    f = f*base;
                }
                f = util::rounding_div_2(f);
            }

            while e > bias
            {
                e = e - U::one();
                f = f*base;
            }
            while e < bias
            {
                e = e + U::one();
                f = util::rounding_div(f, base);
            }

            f
        }
        else
        {
            let base = U::from(EXP_BASE).unwrap();
            let mut f = U::from(f).unwrap();

            for _ in 0..FRAC_SIZE
            {
                while e > bias && f.leading_zeros() > 3
                {
                    e = e - U::one();
                    f = f*base;
                }
                f = util::rounding_div_2(f);
            }

            while e > bias
            {
                e = e - U::one();
                f = f*base;
            }
            while e < bias
            {
                e = e + U::one();
                f = util::rounding_div(f, base);
            }

            match I::from(f)
            {
                Some(f) => f,
                None => I::max_value()
            }
        }
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
        let mut e = self.exp_bits();
        let mut f = self.frac_bits();
        
        if !e.is_zero() || INT_SIZE != 0 //normal
        {
            f = f + (self.int_bits() << Self::INT_POS);
        }
        else
        {
            f = f << 1usize;
        }
        
        let bias = Self::exp_bias();

        let n = if util::bitsize_of::<I>() > util::bitsize_of::<U>()
        {
            let mut f = I::from(f).unwrap();
            let base = I::from(EXP_BASE).unwrap();

            for _ in 0..FRAC_SIZE
            {
                while e > bias && f.leading_zeros() > 3
                {
                    e = e - U::one();
                    f = f*base;
                }
                f = util::rounding_div_2(f);
            }

            while e > bias
            {
                e = e - U::one();
                f = f*base;
            }
            while e < bias
            {
                e = e + U::one();
                f = util::rounding_div(f, base);
            }

            f
        }
        else
        {
            let base = U::from(EXP_BASE).unwrap();
            let mut f = U::from(f).unwrap();

            for _ in 0..FRAC_SIZE
            {
                while e > bias && f.leading_zeros() > 3
                {
                    e = e - U::one();
                    f = f*base;
                }
                f = util::rounding_div_2(f);
            }

            while e > bias
            {
                e = e - U::one();
                f = f*base;
            }
            while e < bias
            {
                e = e + U::one();
                f = util::rounding_div(f, base);
            }

            match I::from(f)
            {
                Some(f) => f,
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
        self = max - ((max - min - ((self - min + Self::one()) % (max - min))) % (max - min)) - Self::one();

        let s = !self.sign_bit().is_zero();
        let mut e = self.exp_bits();
        let mut f = self.frac_bits();
        
        if !e.is_zero() || INT_SIZE != 0 //normal
        {
            f = f + (self.int_bits() << Self::INT_POS);
        }
        else
        {
            f = f << 1usize;
        }
        
        let bias = Self::exp_bias();
        
        let n = if util::bitsize_of::<I>() > util::bitsize_of::<U>()
        {
            let mut f = I::from(f).unwrap();
            let base = I::from(EXP_BASE).unwrap();

            for _ in 0..FRAC_SIZE
            {
                while e > bias && f.leading_zeros() > 3
                {
                    e = e - U::one();
                    f = f*base;
                }
                f = util::rounding_div_2(f);
            }

            while e > bias
            {
                e = e - U::one();
                f = f*base;
            }
            while e < bias
            {
                e = e + U::one();
                f = util::rounding_div(f, base);
            }

            f
        }
        else
        {
            let base = U::from(EXP_BASE).unwrap();
            let mut f = U::from(f).unwrap();

            for _ in 0..FRAC_SIZE
            {
                while e > bias && f.leading_zeros() > 3
                {
                    e = e - U::one();
                    f = f*base;
                }
                f = util::rounding_div_2(f);
            }

            while e > bias
            {
                e = e - U::one();
                f = f*base;
            }
            while e < bias
            {
                e = e + U::one();
                f = util::rounding_div(f, base);
            }

            if s
            {
                if f > U::zero() && f - U::one() >= U::from(-I::one() - I::min_value()).unwrap()
                {
                    return I::min_value()
                }
            }
            else
            {
                if f >= U::from(I::max_value()).unwrap()
                {
                    return I::max_value()
                }
            }
            <I as NumCast>::from(f).unwrap()
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
        (self.to_bits() & U::max_value() >> util::bitsize_of::<U>() - Self::SIGN_POS - Self::SIGN_SIZE) >> Self::SIGN_POS
    }
    pub fn exp_bits(self) -> U
    {
        (self.to_bits() & U::max_value() >> util::bitsize_of::<U>() - Self::EXP_POS - EXP_SIZE) >> Self::EXP_POS
    }
    pub fn int_bits(self) -> U
    {
        if Self::IS_INT_IMPLICIT
        {
            return if self.is_normal() {U::one()} else {U::zero()}
        }
        (self.to_bits() & U::max_value() >> util::bitsize_of::<U>() - Self::INT_POS - INT_SIZE) >> Self::INT_POS
    }
    pub fn frac_bits(self) -> U
    {
        (self.to_bits() & U::max_value() >> util::bitsize_of::<U>() - Self::FRAC_POS - FRAC_SIZE) >> Self::FRAC_POS
    }

    pub fn exp_bias() -> U
    {
        U::max_value() >> util::bitsize_of::<U>() + 1 - EXP_SIZE
    }
    
    pub fn nan() -> Self
    {
        Self::snan()
    }
    
    pub fn qnan() -> Self
    {
        Self::from_bits(U::max_value() >> util::bitsize_of::<U>() - Self::SIGN_POS)
    }
    
    pub fn snan() -> Self
    {
        if Self::INT_POS + INT_SIZE < 1
        {
            return Self::qnan()
        }
        Self::from_bits((U::max_value() >> util::bitsize_of::<U>() - Self::SIGN_POS) - (U::one() << Self::INT_POS + INT_SIZE - 1))
    }

    pub fn is_snan(self) -> bool
    {
        if Self::INT_POS + INT_SIZE < 1 || !self.is_nan()
        {
            return false
        }
        (self.to_bits() >> Self::INT_POS + INT_SIZE - 1) & U::one() == U::zero()
    }

    pub fn infinity() -> Self
    {
        Self::from_bits((U::max_value() >> util::bitsize_of::<U>() - EXP_SIZE) << Self::EXP_POS)
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
        Self::from_bits(((U::max_value() >> util::bitsize_of::<U>() - EXP_SIZE) << Self::EXP_POS) - U::one())
    }

    pub fn is_nan(self) -> bool
    {
        !self.is_finite() && !(self.frac_bits().is_zero() && (INT_SIZE == 0 || self.int_bits().is_zero()))
    }

    pub fn is_infinite(self) -> bool
    {
        !self.is_finite() && self.frac_bits().is_zero() && (INT_SIZE == 0 || self.int_bits().is_zero())
    }

    pub fn is_finite(self) -> bool
    {
        self.exp_bits() != (U::max_value() >> util::bitsize_of::<U>() - EXP_SIZE)
    }

    pub fn is_normal(self) -> bool
    {
        if !self.is_finite()
        {
            return false
        }
        if !Self::IS_INT_IMPLICIT
        {
            return true
        }
        self.exp_bits() != U::zero()
    }

    pub fn classify(self) -> FpCategory
    {
        let e = self.exp_bits();
        if self.is_zero()
        {
            return FpCategory::Zero
        }
        if Self::IS_INT_IMPLICIT && e == U::zero()
        {
            return FpCategory::Subnormal
        }
        if e == (U::max_value() >> util::bitsize_of::<U>() - EXP_SIZE)
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
        if !self.is_finite()
        {
            return self
        }
        let mut m = self % Self::one();
        if self.is_sign_negative() && !m.is_zero()
        {
            m += Self::one()
        }
        self - m
    }

    pub fn ceil(self) -> Self
    {
        if !self.is_finite()
        {
            return self
        }
        let mut m = self % Self::one();
        if self.is_sign_positive() && !m.is_zero()
        {
            m -= Self::one()
        }
        self - m
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
        if !self.is_finite()
        {
            return self
        }
        let m = self % Self::one();
        self - m
    }

    pub fn fract(self) -> Self
    {
        self - self.trunc()
    }

    pub fn abs(self) -> Self
    {
        Self::from_bits(self.to_bits() & (U::max_value() >> util::bitsize_of::<U>() - Self::SIGN_POS))
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

    pub fn powi<I: Int>(self, n: I) -> Self
    {
        util::powi(self, n)
    }
    
    pub fn powu<I: UInt>(self, n: I) -> Self
    {
        util::powu(self, n)
    }

    pub fn powf(self, n: Self) -> Self
    {
        let xabs = self.abs();
        let nabs = n.abs();

        let exp_frac = U::from((1usize << FRAC_SIZE + INT_SIZE).ilog(EXP_BASE)).unwrap();

        let edge_x = {
            let e = (U::one() << EXP_SIZE) - U::one() - U::one();
            let i = if INT_SIZE > 0
            {
                U::one()
            }
            else
            {
                U::zero()
            };
            Self::from_bits((e << Self::EXP_POS) + (i << Self::INT_POS))
        };
        let edge_n = {
            let e = Self::exp_bias() + exp_frac;
            let i = if INT_SIZE > 0
            {
                U::one()
            }
            else
            {
                U::zero()
            };
            Self::from_bits((e << Self::EXP_POS) + (i << Self::INT_POS))
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
        let n_d = (nabs + Self::one())*<Self as From<_>>::from(0.5);
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

        if (xabs >= edge_x) || (nabs >= edge_n)
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
                let i: u128 = nabs.to_uint().unwrap();
                
                let x = if s
                {
                    Self::one() / self
                }
                else
                {
                    self
                };
                
                return x.powu(i)
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

        let xabs_log = xabs.log_base();

        let n_xabs_log = n*xabs_log;
        
        n_xabs_log.exp_base()
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

        let y = if EXP_BASE != 2
        {
            let xabs_log = self.log_base();

            let n_xabs_log = <Self as From<_>>::from(0.5)*xabs_log;
            
            n_xabs_log.exp_base()
        }
        else if !Self::IS_INT_IMPLICIT
        {
            let x = Fp::<U, EXP_SIZE, 0, FRAC_SIZE, EXP_BASE>::from_fp(self);
            Self::from_fp(Fp::<U, EXP_SIZE, 0, FRAC_SIZE, EXP_BASE>::from_bits(
                ((Fp::<U, EXP_SIZE, 0, FRAC_SIZE, EXP_BASE>::exp_bias() + U::one()) * (U::one() << (FRAC_SIZE - 1)))
                + (x.to_bits() >> 1usize)
                - (U::one() << (FRAC_SIZE - 1))
            ))
        }
        else
        {
            Self::from_bits(
                ((Self::exp_bias() + U::one()) * (U::one() << (FRAC_SIZE - 1)))
                + (self.to_bits() >> 1usize)
                - (U::one() << (FRAC_SIZE - 1))
            )
        };

        const NEWTON: usize = 4;
        let half = <Self as From<_>>::from(0.5);
        let mut y = y;
        for _ in 0..NEWTON
        {
            y = half*(y + self/y);
        }
        y
    }

    pub fn exp_base(self) -> Self
    {
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
        
        let exp_frac = U::from((1usize << FRAC_SIZE + INT_SIZE).ilog(EXP_BASE) as usize).unwrap();

        if self >= Self::from_uint(U::one() << (EXP_SIZE - 1))
        {
            return Self::from_bits(((U::one() << EXP_SIZE) - U::one()) << Self::EXP_POS)*self;
        }
        if self <= -Self::from_uint((U::one() << (EXP_SIZE - 1)) + exp_frac - U::one())
        {
            return -Self::from_bits(U::one() << (Self::EXP_POS - EXP_SIZE/2))*self;
        }

        let neg = self.is_sign_negative();
        let x = self.abs().max(-Self::from_uint(Self::exp_bias() - U::one()));
        
        let e = x.floor();
        let f = x - e;

        let z = <Self as From<_>>::from((Into::<f64>::into(f)*(EXP_BASE as f64).ln()).exp());

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
        let i = if INT_SIZE > 0
        {
            U::one()
        }
        else
        {
            U::zero()
        };
        let y = Self::from_bits((e << Self::EXP_POS) + (i << Self::INT_POS))*z;

        if neg
        {
            return y.recip()
        }
        y
    }

    pub fn exp(self) -> Self
    {
        if EXP_BASE == 2
        {
            return (self/Self::LN_2()).exp_base()
        }
        if EXP_BASE == 10
        {
            return (self/Self::LN_10()).exp_base()
        }
        (self/<Self as From<_>>::from((EXP_BASE as f64).ln())).exp_base()
    }

    pub fn exp10(self) -> Self
    {
        if EXP_BASE == 10
        {
            return self.exp_base()
        }
        (self*Self::LN_10()).exp()
    }

    pub fn exp2(self) -> Self
    {
        if EXP_BASE == 2
        {
            return self.exp_base()
        }
        (self*Self::LN_2()).exp()
    }

    pub fn ln(self) -> Self
    {
        if EXP_BASE == 2
        {
            return self.log_base()*Self::LN_2()
        }
        if EXP_BASE == 10
        {
            return self.log_base()*Self::LN_10()
        }
        self.log_base()*<Self as From<_>>::from((EXP_BASE as f64).ln())
    }

    pub fn log(self, base: Self) -> Self
    {
        self.ln()/base.ln()
    }

    pub fn log_base(self) -> Self
    {
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

        let e = self.exp_bits();
        let bias = Self::exp_bias();
        let mut y = Self::from_uint(e)
            - Self::from_uint(bias);

        let mut f = self.frac_bits();
        if !e.is_zero() || !Self::IS_INT_IMPLICIT //normal
        {
            f = f + (self.int_bits() << Self::INT_POS);
        }
        else
        {
            f = f << 1usize;
        }

        if Self::IS_INT_IMPLICIT
        {
            f = f - (U::one() << FRAC_SIZE);
        }

        let u = Self::from_bits(f + (bias << Self::EXP_POS));
        
        let u: f64 = u.into();
        y += <Self as From<_>>::from(u.log(EXP_BASE as f64)); 
        y
    }

    pub fn log2(self) -> Self
    {
        if EXP_BASE == 2
        {
            return self.log_base()
        }
        self.ln()/Self::LN_2()
    }

    pub fn log10(self) -> Self
    {
        if EXP_BASE == 10
        {
            return self.log_base()
        }
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
        let bias = Self::exp_bias();
        
        if !Self::IS_INT_IMPLICIT
        {
            return Self::from_bits((bias << Self::EXP_POS) + U::one())
        }

        let exp_frac = U::from((1usize << FRAC_SIZE).ilog(EXP_BASE) as usize).unwrap();
        if bias <= exp_frac
        {
            return Self::from_bits(util::powu(U::from(EXP_BASE).unwrap(), <usize as NumCast>::from(bias).unwrap() - 1))
        }
        Self::from_bits((bias - exp_frac) << Self::EXP_POS)
    }

    pub fn copysign(self, sign: Self) -> Self
    {
        let s = sign.sign_bit();
        let e = self.exp_bits();
        let f = self.frac_bits();
        if !Self::IS_INT_IMPLICIT
        {
            let i = self.int_bits();
            return Self::from_bits((s << Self::SIGN_POS) + (e << Self::EXP_POS) + (i << Self::INT_POS) + (f << Self::FRAC_POS))
        }
        Self::from_bits((s << Self::SIGN_POS) + (e << Self::EXP_POS) + (f << Self::FRAC_POS))
    }

    pub fn is_subnormal(self) -> bool
    {
        Self::IS_INT_IMPLICIT && self.exp_bits() == U::zero() && !self.is_zero()
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
        if !self.is_finite()
        {
            return false;
        }
        if !Self::IS_INT_IMPLICIT
        {
            return self.frac_bits().is_zero() && self.int_bits().is_zero()
        }
        self.exp_bits().is_zero() && self.frac_bits().is_zero()
    }

    fn set_zero(&mut self)
    {
        *self = Self::zero()
    }
    
    pub fn one() -> Self
    {
        let bias = Self::exp_bias();
        if !Self::IS_INT_IMPLICIT
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
        let s0 = self.sign_bit();
        let s1 = other.sign_bit();

        if s0 != s1
        {
            return s1.cmp(&s0)
        }

        let s = !s0.is_zero();
        
        let c0 = if self.is_nan()
        {
            if !self.is_snan()
            {
                3
            }
            else
            {
                2
            }
        }
        else if !self.is_finite()
        {
            1
        }
        else
        {
            0
        };
        
        let c1 = if other.is_nan()
        {
            if !other.is_snan()
            {
                3
            }
            else
            {
                2
            }
        }
        else if !other.is_finite()
        {
            1
        }
        else
        {
            0
        };

        if c0 != 0 || c1 != 0
        {
            return if s {c1.cmp(&c0)} else {c0.cmp(&c1)}
        }

        let mut e0 = self.exp_bits();
        let mut e1 = other.exp_bits();

        let mut f0 = self.frac_bits();
        let mut f1 = other.frac_bits();

        if !Self::IS_INT_IMPLICIT
        {
            let base = U::from(EXP_BASE).unwrap();

            f0 = f0 + (self.int_bits() << FRAC_SIZE);
            f1 = f1 + (other.int_bits() << FRAC_SIZE);
            
            while e0 > e1 && f0 < U::one() << Self::MANTISSA_OP_SIZE - Self::BASE_SIZE
            {
                e0 = e0 - U::one();
                f0 = f0*base;
            }            
            while e1 > e0 && f1 < U::one() << Self::MANTISSA_OP_SIZE - Self::BASE_SIZE
            {
                e1 = e1 - U::one();
                f1 = f1*base;
            }
            
            while e0 < e1 && f0 >= U::one() << Self::MANTISSA_OP_SIZE
            {
                e0 = e0 + U::one();
                f0 = util::rounding_div(f0, base);
            }
            while e1 < e0 && f1 >= U::one() << Self::MANTISSA_OP_SIZE
            {
                e1 = e1 + U::one();
                f1 = util::rounding_div(f1, base);
            }
        }
        
        if e0 != e1
        {
            return if s {e1.cmp(&e0)} else {e0.cmp(&e1)}
        }

        if s {f1.cmp(&f0)} else {f0.cmp(&f1)}
    }

    pub fn mul_base(self) -> Self
    {
        if !self.is_finite()
        {
            return self
        }
        let s = self.sign_bit();
        let mut e = self.exp_bits();
        let mut f = self.frac_bits();
        
        if !e.is_zero() || !Self::IS_INT_IMPLICIT //normal
        {
            f = f + (self.int_bits() << FRAC_SIZE);
        }
        else
        {
            f = f << 1usize
        }

        e = e + U::one();

        let base = U::from(EXP_BASE).unwrap();        
        while e > U::zero() && f < U::one() << Self::MANTISSA_OP_SIZE - Self::BASE_SIZE
        {
            e = e - U::one();
            f = f*base;
        }
        while e < U::one() << EXP_SIZE && f >= U::one() << Self::MANTISSA_OP_SIZE
        {
            e = e + U::one();
            f = util::rounding_div(f, base);
        }

        if e >= (U::one() << EXP_SIZE) - U::one()
        {
            return if !s.is_zero() {Self::neg_infinity()} else {Self::infinity()}
        }

        let s_bit = s << Self::SIGN_POS;
        if e.is_zero() && Self::IS_INT_IMPLICIT // subnormal
        {
            Fp::from_bits(s_bit + (f >> 1usize))
        }
        else
        {
            if Self::IS_INT_IMPLICIT
            {
                f = f - (U::one() << FRAC_SIZE);
            }

            if e >= (U::one() << EXP_SIZE) - U::one()
            {
                return if !s.is_zero() {Self::neg_infinity()} else {Self::infinity()}
            }

            Fp::from_bits(s_bit + f + (e << Self::EXP_POS))
        }
    }

    pub fn div_base(self) -> Self
    {
        if !self.is_finite()
        {
            return self
        }
        let s = self.sign_bit();
        let mut e = self.exp_bits();
        let mut f = self.frac_bits();
        
        if !e.is_zero() || !Self::IS_INT_IMPLICIT //normal
        {
            f = f + (self.int_bits() << FRAC_SIZE);
        }
        else
        {
            f = f << 1usize
        }

        let base = U::from(EXP_BASE).unwrap();  
        if e > U::zero()
        {
            e = e - U::one();
        }
        else
        {
            f = util::rounding_div(f, base);
        }
              
        while e > U::zero() && f < U::one() << Self::MANTISSA_OP_SIZE - Self::BASE_SIZE
        {
            e = e - U::one();
            f = f*base;
        }
        while e < U::one() << EXP_SIZE && f >= U::one() << Self::MANTISSA_OP_SIZE
        {
            e = e + U::one();
            f = util::rounding_div(f, base);
        }

        if e >= (U::one() << EXP_SIZE) - U::one()
        {
            return if !s.is_zero() {Self::neg_infinity()} else {Self::infinity()}
        }

        let s_bit = s << Self::SIGN_POS;
        if e.is_zero() && Self::IS_INT_IMPLICIT // subnormal
        {
            Fp::from_bits(s_bit + (f >> 1usize))
        }
        else
        {
            if Self::IS_INT_IMPLICIT
            {
                f = f - (U::one() << FRAC_SIZE);
            }

            if e >= (U::one() << EXP_SIZE) - U::one()
            {
                return if !s.is_zero() {Self::neg_infinity()} else {Self::infinity()}
            }

            Fp::from_bits(s_bit + f + (e << Self::EXP_POS))
        }
    }
}

impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:,
    [(); EXP_BASE - 2]:,
    U: ConstZero
{
    pub const ZERO: Self = Self::from_bits(U::ZERO);
}

impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:,
    [(); EXP_BASE - 2]:,
    U: ConstOne
{
    pub const MIN_POSITIVE: Self = Self::from_bits(U::ONE);
}