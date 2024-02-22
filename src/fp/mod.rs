use std::cmp::Ordering;
use std::num::FpCategory;

use array_math::{ArrayMath, ArrayOps};
use num_traits::{ConstOne, ConstZero, FloatConst, FromBytes, ToBytes};
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
        from,
        into,
        num,
        partial_eq,
        partial_ord,
        signed,
    }
);

const NEWTON_EXP: usize = 5;
const NEWTON_LN: usize = 4;
const NEWTON_RT: usize = 4;
const NEWTON_TRIG: usize = 4;

/// A custom floating point type.
/// 
/// The bit layout is as follows:
/// ```txt
/// No data: | Sign: | Exponent:  | Integer:   | Fractional: |
/// <  ..  > | < 1 > | <EXP_SIZE> | <INT_SIZE> | <FRAC_SIZE> |
/// ```
/// 
/// The value of a real floating-point number is the following:
/// ```txt
/// x = (-1)**sign*EXP_BASE**(exponent - bias)*mantissa
/// ```
/// 
/// where the bias equals
/// ```txt
/// bias = 2**(EXP_SIZE - 1) - 1
/// ```
///
/// If the exponent has the maximum value, the number is either infinity or NaN.
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
    /// Size of floating-point number in bits
    pub const BIT_SIZE: usize = EXP_SIZE + INT_SIZE + FRAC_SIZE + 1;
    /// Size of the sign bit
    pub const SIGN_SIZE: usize = 1;

    /// Position of the sign bit
    pub const SIGN_POS: usize = EXP_SIZE + INT_SIZE + FRAC_SIZE;
    /// Position of the first exponent bit
    pub const EXP_POS: usize = INT_SIZE + FRAC_SIZE;
    /// Position of the first integer bit
    pub const INT_POS: usize = FRAC_SIZE;
    /// Position of the first fractional bit
    pub const FRAC_POS: usize = 0;

    /// Number of significant digits in base 2.
    pub const MANTISSA_DIGITS: usize = INT_SIZE + FRAC_SIZE;

    /// `true` if the number contains an implicit integer bit
    pub const IS_INT_IMPLICIT: bool = INT_SIZE == 0;

    const MANTISSA_OP_SIZE: usize = FRAC_SIZE + INT_SIZE + Self::IS_INT_IMPLICIT as usize;
    const BASE_PADDING: usize = util::bitsize_of::<usize>() - EXP_BASE.leading_zeros() as usize - 1;

    /// Converts from one custom floating-point number to another.
    /// Rounding errors may occurr.
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn from_fp<V: UInt, const E: usize, const I: usize, const F: usize, const B: usize>(fp: Fp<V, E, I, F, B>) -> Self
    where
        [(); util::bitsize_of::<V>() - E - I - F - 1]:,
        [(); util::bitsize_of::<V>() - E - 0 - F - 1]:,
        [(); B - 2]:
    {
        if EXP_SIZE == E && INT_SIZE == I && EXP_BASE == B
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
                        util::rounding_div_pow(b, U::from(2).unwrap(), F - FRAC_SIZE)
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
                    util::rounding_div_pow(fp.to_bits(), V::from(2).unwrap(), F - FRAC_SIZE)
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
            while f.leading_zeros() as usize <= Fp::<V, E, I, F, B>::BASE_PADDING
            {
                e = e + U::one();
                f = util::rounding_div(f, base2);
            }
            f = f*base1
        }
        while e1 < bias1
        {
            e1 = e1 + V::one();
            while e > U::zero() && f.leading_zeros() as usize > Self::BASE_PADDING
            {
                e = e - U::one();
                f = f*base2;
            }
            f = util::rounding_div(f, base1)
        }
        
        let s_bit = s << Self::SIGN_POS;

        while e > U::zero() && f < U::one() << Self::MANTISSA_OP_SIZE - Self::BASE_PADDING
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

    /// Converts an unsigned integer into a custom floating-point type.
    #[must_use = "method returns a new number and does not mutate the original value"]
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

            while e > U::zero() && f < I::one() << Self::MANTISSA_OP_SIZE - Self::BASE_PADDING
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
            
            while e > U::zero() && f < U::one() << Self::MANTISSA_OP_SIZE - Self::BASE_PADDING
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
    
    /// Converts a signed integer into a custom floating-point type.
    #[must_use = "method returns a new number and does not mutate the original value"]
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

            while e > U::zero() && f < I::one() << Self::MANTISSA_OP_SIZE - Self::BASE_PADDING
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

            while e > U::zero() && f < U::one() << Self::MANTISSA_OP_SIZE - Self::BASE_PADDING
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

    /// Converts a custom floating-point type into an unsigned integer.
    ///
    /// Returns [None](None) if out of bounds.
    #[must_use = "method returns a new number and does not mutate the original value"]
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
    
    /// Converts a custom floating-point type into an unsigned integer.
    ///
    /// Wraps if out of bounds.
    #[must_use = "method returns a new number and does not mutate the original value"]
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
    
    /// Converts a custom floating-point type into a signed integer.
    ///
    /// Returns [None](None) if out of bounds.
    #[must_use = "method returns a new number and does not mutate the original value"]
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
    
    /// Converts a custom floating-point type into a signed integer.
    ///
    /// Wraps if out of bounds.
    #[must_use = "method returns a new number and does not mutate the original value"]
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

    /// Raw transmutation from bits.
    ///
    /// Note that this function is distinct from [`Fp::from_uint`], which attempts to
    /// preserve the *numeric* value, and not the bitwise value.
    ///
    /// # Examples
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    ///
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let v = FpSingle::from_bits(0x41480000);
    /// assert_eq!(v, FpSingle::from(12.5));
    /// ```
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub const fn from_bits(bits: U) -> Self
    {
        Self(bits)
    }

    /// Raw transmutation to bits.
    ///
    /// Note that this function is distinct from [`Fp::to_uint`], which attempts to
    /// preserve the *numeric* value, and not the bitwise value.
    ///
    /// # Examples
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    ///
    /// use custom_float::ieee754::FpSingle;
    ///
    /// assert_ne!(FpSingle::from(1.0).to_bits(), FpSingle::from(1.0).to_uint().unwrap()); // to_bits() is not casting!
    /// assert_eq!(FpSingle::from(12.5).to_bits(), 0x41480000);
    /// ```
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub const fn to_bits(self) -> U
    {
        self.0
    }
    
    /// Return the memory representation of this floating point number as a byte array in
    /// big-endian (network) byte order.
    ///
    /// # Examples
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    ///
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let bytes = FpSingle::from(12.5).to_be_bytes();
    /// assert_eq!(bytes, [0x41, 0x48, 0x00, 0x00]);
    /// ```
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub fn to_be_bytes(self) -> U::Bytes
    where
        U: ToBytes
    {
        self.to_bits().to_be_bytes()
    }
    
    /// Return the memory representation of this floating point number as a byte array in
    /// little-endian byte order.
    ///
    /// # Examples
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    ///
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let bytes = FpSingle::from(12.5).to_le_bytes();
    /// assert_eq!(bytes, [0x00, 0x00, 0x48, 0x41]);
    /// ```
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub fn to_le_bytes(self) -> U::Bytes
    where
        U: ToBytes
    {
        self.to_bits().to_le_bytes()
    }
    
    /// Return the memory representation of this floating point number as a byte array in
    /// native byte order.
    ///
    /// As the target platform's native endianness is used, portable code
    /// should use [`to_be_bytes`] or [`to_le_bytes`], as appropriate, instead.
    ///
    /// [`to_be_bytes`]: Self::to_be_bytes
    /// [`to_le_bytes`]: Self::to_le_bytes
    ///
    /// # Examples
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    ///
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let bytes = FpSingle::from(12.5).to_ne_bytes();
    /// assert_eq!(
    ///     bytes,
    ///     if cfg!(target_endian = "big") {
    ///         [0x41, 0x48, 0x00, 0x00]
    ///     } else {
    ///         [0x00, 0x00, 0x48, 0x41]
    ///     }
    /// );
    /// ```
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub fn to_ne_bytes(self) -> U::Bytes
    where
        U: ToBytes
    {
        self.to_bits().to_ne_bytes()
    }
    
    /// Create a floating point value from its representation as a byte array in big endian.
    ///
    /// # Examples
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    ///
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let value = FpSingle::from_be_bytes(&[0x41, 0x48, 0x00, 0x00]);
    /// assert_eq!(value, FpSingle::from(12.5));
    /// ```
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub fn from_be_bytes(bytes: &U::Bytes) -> Self
    where
        U: FromBytes
    {
        Self::from_bits(U::from_be_bytes(bytes))
    }

    /// Create a floating point value from its representation as a byte array in little endian.
    ///
    /// See [`from_bits`](Self::from_bits) for some discussion of the
    /// portability of this operation (there are almost no issues).
    ///
    /// # Examples
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    ///
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let value = FpSingle::from_le_bytes(&[0x00, 0x00, 0x48, 0x41]);
    /// assert_eq!(value, FpSingle::from(12.5));
    /// ```
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub fn from_le_bytes(bytes: &U::Bytes) -> Self
    where
        U: FromBytes
    {
        Self::from_bits(U::from_le_bytes(bytes))
    }

    /// Create a floating point value from its representation as a byte array in native endian.
    ///
    /// As the target platform's native endianness is used, portable code
    /// likely wants to use [`from_be_bytes`] or [`from_le_bytes`], as
    /// appropriate instead.
    ///
    /// [`from_be_bytes`]: f32::from_be_bytes
    /// [`from_le_bytes`]: f32::from_le_bytes
    ///
    /// See [`from_bits`](Self::from_bits) for some discussion of the
    /// portability of this operation (there are almost no issues).
    ///
    /// # Examples
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    ///
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let value = FpSingle::from_ne_bytes(if cfg!(target_endian = "big") {
    ///     &[0x41, 0x48, 0x00, 0x00]
    /// } else {
    ///     &[0x00, 0x00, 0x48, 0x41]
    /// });
    /// assert_eq!(value, FpSingle::from(12.5));
    /// ```
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub fn from_ne_bytes(bytes: &U::Bytes) -> Self
    where
        U: FromBytes
    {
        Self::from_bits(U::from_ne_bytes(bytes))
    }

    /// Returns the sign bit of the custom floating-point number.
    #[must_use = "this returns the result of the operation, without modifying the original"]
    pub fn sign_bit(self) -> U
    {
        (self.to_bits() & U::max_value() >> util::bitsize_of::<U>() - Self::SIGN_POS - Self::SIGN_SIZE) >> Self::SIGN_POS
    }
    /// Returns the exponent bits of the custom floating-point number.
    #[must_use = "this returns the result of the operation, without modifying the original"]
    pub fn exp_bits(self) -> U
    {
        (self.to_bits() & U::max_value() >> util::bitsize_of::<U>() - Self::EXP_POS - EXP_SIZE) >> Self::EXP_POS
    }
    /// Returns the integer bits of the custom floating-point number.
    #[must_use = "this returns the result of the operation, without modifying the original"]
    pub fn int_bits(self) -> U
    {
        if Self::IS_INT_IMPLICIT
        {
            return if self.is_normal() {U::one()} else {U::zero()}
        }
        (self.to_bits() & U::max_value() >> util::bitsize_of::<U>() - Self::INT_POS - INT_SIZE) >> Self::INT_POS
    }
    /// Returns the fractional bits of the custom floating-point number.
    #[must_use = "this returns the result of the operation, without modifying the original"]
    pub fn frac_bits(self) -> U
    {
        (self.to_bits() & U::max_value() >> util::bitsize_of::<U>() - Self::FRAC_POS - FRAC_SIZE) >> Self::FRAC_POS
    }

    /// Returns the exponent bias
    #[must_use]
    pub fn exp_bias() -> U
    {
        U::max_value() >> util::bitsize_of::<U>() + 1 - EXP_SIZE
    }
    
    /// Returns the `NaN` value.
    ///
    /// ```
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let nan = FpSingle::nan();
    ///
    /// assert!(nan.is_nan());
    /// ```
    #[must_use]
    pub fn nan() -> Self
    {
        Self::snan()
    }
    
    /// Returns the `qNaN` value.
    ///
    /// ```
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let qnan = FpSingle::qnan();
    ///
    /// assert!(qnan.is_nan());
    /// assert!(!qnan.is_snan());
    /// ```
    #[must_use]
    pub fn qnan() -> Self
    {
        Self::from_bits(U::max_value() >> util::bitsize_of::<U>() - Self::SIGN_POS)
    }
    
    /// Returns the `sNaN` value.
    ///
    /// ```
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let snan = FpSingle::snan();
    ///
    /// assert!(snan.is_nan());
    /// assert!(snan.is_snan());
    /// ```
    #[must_use]
    pub fn snan() -> Self
    {
        if Self::INT_POS + INT_SIZE < 1
        {
            return Self::qnan()
        }
        Self::from_bits((U::max_value() >> util::bitsize_of::<U>() - Self::SIGN_POS) - (U::one() << Self::INT_POS + INT_SIZE - 1))
    }

    /// Returns `true` if the number is a signaling NaN.
    ///
    /// ```
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let snan = FpSingle::snan();
    /// let qnan = FpSingle::qnan();
    ///
    /// assert!(snan.is_snan());
    /// assert!(!qnan.is_snan());
    /// ```
    #[must_use = "this returns the result of the operation, without modifying the original"]
    pub fn is_snan(self) -> bool
    {
        if Self::INT_POS + INT_SIZE < 1 || !self.is_nan()
        {
            return false
        }
        (self.to_bits() >> Self::INT_POS + INT_SIZE - 1) & U::one() == U::zero()
    }

    /// Returns the infinite value.
    ///
    /// ```
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let infinity = FpSingle::infinity();
    ///
    /// assert!(infinity.is_infinite());
    /// assert!(!infinity.is_finite());
    /// assert!(infinity > FpSingle::max_value());
    /// ```
    #[must_use]
    pub fn infinity() -> Self
    {
        Self::from_bits((U::max_value() >> util::bitsize_of::<U>() - EXP_SIZE) << Self::EXP_POS)
    }

    /// Returns the negative infinite value.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    ///
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let neg_infinity = FpSingle::neg_infinity();
    ///
    /// assert!(neg_infinity.is_infinite());
    /// assert!(!neg_infinity.is_finite());
    /// assert!(neg_infinity < FpSingle::min_value());
    /// ```
    #[must_use]
    pub fn neg_infinity() -> Self
    {
        -Self::infinity()
    }

    /// Returns `-0.0`.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let inf = FpSingle::infinity();
    /// let zero = FpSingle::zero();
    /// let neg_zero = FpSingle::neg_zero();
    ///
    /// assert_eq!(zero, neg_zero);
    /// assert_eq!(FpSingle::from(7.0)/inf, zero);
    /// assert_eq!(zero * FpSingle::from(10.0), zero);
    /// ```
    #[must_use]
    pub fn neg_zero() -> Self
    {
        -Self::zero()
    }

    /// Returns the smallest finite value.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use std::f64;
    ///
    /// let x = FpDouble::min_value();
    ///
    /// assert_eq!(x, FpDouble::from(f64::MIN));
    /// ```
    #[must_use]
    pub fn min_value() -> Self
    {
        -Self::max_value()
    }

    /// Returns the smallest positive, normal value.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use std::f64;
    ///
    /// let x = FpDouble::min_positive_value();
    ///
    /// assert_eq!(x, FpDouble::from(f64::MIN_POSITIVE));
    /// ```
    #[must_use]
    pub fn min_positive_value() -> Self
    {
        if Self::IS_INT_IMPLICIT
        {
            Self::from_bits(U::one() << Self::EXP_POS)
        }
        else
        {
            Self::from_bits((U::one() << Self::EXP_POS) + (U::one() << Self::INT_POS))
        }
    }

    /// [Machine epsilon] value.
    ///
    /// This is the difference between `1.0` and the next larger representable number.
    ///
    /// [Machine epsilon]: https://en.wikipedia.org/wiki/Machine_epsilon
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use std::f64;
    ///
    /// let x = FpDouble::epsilon();
    ///
    /// assert_eq!(x, FpDouble::from(f64::EPSILON));
    /// ```
    #[must_use]
    pub fn epsilon() -> Self
    {
        let bias = Self::exp_bias();
        
        if !Self::IS_INT_IMPLICIT
        {
            return Self::from_bits((bias << Self::EXP_POS) + U::one())
        }

        let exp_frac = U::from(util::count_digits_in_base(FRAC_SIZE + INT_SIZE, EXP_BASE)).unwrap();
        if bias <= exp_frac
        {
            return Self::from_bits(util::powu(U::from(EXP_BASE).unwrap(), <usize as NumCast>::from(bias).unwrap() - 1))
        }
        Self::from_bits((bias - exp_frac) << Self::EXP_POS)
    }

    /// Returns the largest finite value.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use std::f64;
    ///
    /// let x = FpDouble::max_value();
    /// assert_eq!(x, FpDouble::from(f64::MAX));
    /// ```
    #[must_use]
    pub fn max_value() -> Self
    {
        Self::from_bits(((U::max_value() >> util::bitsize_of::<U>() - EXP_SIZE) << Self::EXP_POS) - U::one())
    }

    /// Returns `true` if this value is NaN.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let nan = FpDouble::nan();
    /// let f = FpDouble::from(7.0);
    ///
    /// assert!(nan.is_nan());
    /// assert!(!f.is_nan());
    /// ```
    #[must_use = "this returns the result of the operation, without modifying the original"]
    pub fn is_nan(self) -> bool
    {
        !self.is_finite() && !(self.frac_bits().is_zero() && (INT_SIZE == 0 || self.int_bits().is_zero()))
    }

    /// Returns `true` if this value is positive infinity or negative infinity, and
    /// `false` otherwise.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let f = FpSingle::from(7.0);
    /// let inf = FpSingle::infinity();
    /// let neg_inf = FpSingle::neg_infinity();
    /// let nan = FpSingle::nan();
    ///
    /// assert!(!f.is_infinite());
    /// assert!(!nan.is_infinite());
    ///
    /// assert!(inf.is_infinite());
    /// assert!(neg_inf.is_infinite());
    /// ```
    #[must_use = "this returns the result of the operation, without modifying the original"]
    pub fn is_infinite(self) -> bool
    {
        !self.is_finite() && self.frac_bits().is_zero() && (INT_SIZE == 0 || self.int_bits().is_zero())
    }

    /// Returns `true` if this number is neither infinite nor NaN.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let f = FpSingle::from(7.0);
    /// let inf = FpSingle::infinity();
    /// let neg_inf = FpSingle::neg_infinity();
    /// let nan = FpSingle::nan();
    ///
    /// assert!(f.is_finite());
    ///
    /// assert!(!nan.is_finite());
    /// assert!(!inf.is_finite());
    /// assert!(!neg_inf.is_finite());
    /// ```
    #[must_use = "this returns the result of the operation, without modifying the original"]
    pub fn is_finite(self) -> bool
    {
        self.exp_bits() != (U::max_value() >> util::bitsize_of::<U>() - EXP_SIZE)
    }

    /// Returns `true` if the number is neither zero, infinite,
    /// [subnormal], or NaN.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let min = FpSingle::min_positive_value(); // 1.17549435e-38f32
    /// let max = FpSingle::max_value();
    /// let lower_than_min = FpSingle::from(1.0e-40_f32);
    /// let zero = FpSingle::zero();
    ///
    /// assert!(min.is_normal());
    /// assert!(max.is_normal());
    ///
    /// assert!(!zero.is_normal());
    /// assert!(!FpSingle::nan().is_normal());
    /// assert!(!FpSingle::infinity().is_normal());
    /// // Values between `0` and `min` are Subnormal.
    /// assert!(!lower_than_min.is_normal());
    /// ```
    /// [subnormal]: http://en.wikipedia.org/wiki/Denormal_number
    #[must_use = "this returns the result of the operation, without modifying the original"]
    pub fn is_normal(self) -> bool
    {
        if !self.is_finite()
        {
            return false
        }
        if !Self::IS_INT_IMPLICIT
        {
            return !self.is_zero()
        }
        self.exp_bits() != U::zero()
    }

    /// Returns `true` if the number is [subnormal].
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let min = FpDouble::min_positive_value(); // 2.2250738585072014e-308_f64
    /// let max = FpDouble::max_value();
    /// let lower_than_min = FpDouble::from(1.0e-308_f64);
    /// let zero = FpDouble::zero();
    ///
    /// assert!(!min.is_subnormal());
    /// assert!(!max.is_subnormal());
    ///
    /// assert!(!zero.is_subnormal());
    /// assert!(!FpDouble::nan().is_subnormal());
    /// assert!(!FpDouble::infinity().is_subnormal());
    /// // Values between `0` and `min` are Subnormal.
    /// assert!(lower_than_min.is_subnormal());
    /// ```
    /// [subnormal]: https://en.wikipedia.org/wiki/Denormal_number
    #[must_use = "this returns the result of the operation, without modifying the original"]
    pub fn is_subnormal(self) -> bool
    {
        Self::IS_INT_IMPLICIT && self.exp_bits() == U::zero() && !self.is_zero()
    }

    /// Returns the floating point category of the number. If only one property
    /// is going to be tested, it is generally faster to use the specific
    /// predicate instead.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use std::num::FpCategory;
    ///
    /// let num = FpDouble::from(12.4f32);
    /// let inf = FpDouble::infinity();
    ///
    /// assert_eq!(num.classify(), FpCategory::Normal);
    /// assert_eq!(inf.classify(), FpCategory::Infinite);
    /// ```
    #[must_use = "this returns the result of the operation, without modifying the original"]
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

    /// Returns the largest integer less than or equal to `self`.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let f = FpDouble::from(3.99);
    /// let g = FpDouble::from(3.0);
    ///
    /// assert_eq!(f.floor(), FpDouble::from(3.0));
    /// assert_eq!(g.floor(), FpDouble::from(3.0));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
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

    /// Returns the smallest integer greater than or equal to `self`.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let f = FpDouble::from(3.01);
    /// let g = FpDouble::from(4.0);
    ///
    /// assert_eq!(f.ceil(), FpDouble::from(4.0));
    /// assert_eq!(g.ceil(), FpDouble::from(4.0));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
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

    /// Returns the nearest integer to `self`. If a value is half-way between two
    /// integers, round away from `0.0`.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let f = FpDouble::from(3.3);
    /// let g = FpDouble::from(-3.3);
    ///
    /// assert_eq!(f.round(), FpDouble::from(3.0));
    /// assert_eq!(g.round(), FpDouble::from(-3.0));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn round(self) -> Self
    {
        if !self.is_finite()
        {
            return self
        }
        let mut m = self % Self::one();
        let half = <Self as From<_>>::from(0.5);
        if self.is_sign_positive() && m >= half
        {
            m -= Self::one()
        }
        if self.is_sign_negative() && m <= -half
        {
            m += Self::one()
        }
        self - m
    }
    
    /// Returns the nearest integer to a number. Rounds half-way cases to the number
    /// with an even least significant digit.
    ///
    /// # Examples
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let f = FpSingle::from(3.3);
    /// let g = FpSingle::from(-3.3);
    /// let h = FpSingle::from(3.5);
    /// let i = FpSingle::from(4.5);
    ///
    /// assert_eq!(f.round_ties_even(), FpSingle::from(3.0));
    /// assert_eq!(g.round_ties_even(), FpSingle::from(-3.0));
    /// assert_eq!(h.round_ties_even(), FpSingle::from(4.0));
    /// assert_eq!(i.round_ties_even(), FpSingle::from(4.0));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn round_ties_even(self) -> Self
    {
        if !self.is_finite()
        {
            return self
        }
        let one = Self::one();
        let t = self % <Self as From<_>>::from(2.0);
        let mut m = self % one;
        let half = <Self as From<_>>::from(0.5);
        if self.is_sign_positive() && m >= half
        {
            if m != half || t > one
            {
                m -= Self::one()
            }
        }
        if self.is_sign_negative() && m <= -half
        {
            if m != -half || t < -one
            {
                m += Self::one()
            }
        }
        self - m
    }

    /// Returns the integer part of `self`.
    /// This means that non-integer numbers are always truncated towards zero.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let f = FpDouble::from(3.3);
    /// let g = FpDouble::from(-3.7);
    ///
    /// assert_eq!(f.trunc(), FpDouble::from(3.0));
    /// assert_eq!(g.trunc(), FpDouble::from(-3.0));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn trunc(self) -> Self
    {
        if !self.is_finite()
        {
            return self
        }
        let m = self % Self::one();
        self - m
    }

    /// Returns the fractional part of a number.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let x = FpDouble::from(3.5);
    /// let y = FpDouble::from(-3.5);
    /// let abs_difference_x = (x.fract() - FpDouble::from(0.5)).abs();
    /// let abs_difference_y = (y.fract() - FpDouble::from(-0.5)).abs();
    ///
    /// assert!(abs_difference_x < FpDouble::from(1e-10));
    /// assert!(abs_difference_y < FpDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn fract(self) -> Self
    {
        self - self.trunc()
    }

    /// Computes the absolute value of `self`. Returns `NaN` if the number is `NaN`.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let x = FpDouble::from(3.5);
    /// let y = FpDouble::from(-3.5);
    ///
    /// let abs_difference_x = (x.abs() - x).abs();
    /// let abs_difference_y = (y.abs() - (-y)).abs();
    ///
    /// assert!(abs_difference_x < FpDouble::from(1e-10));
    /// assert!(abs_difference_y < FpDouble::from(1e-10));
    ///
    /// assert!(FpDouble::nan().abs().is_nan());
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn abs(self) -> Self
    {
        Self::from_bits(self.to_bits() & (U::max_value() >> util::bitsize_of::<U>() - Self::SIGN_POS))
    }

    /// Returns a number that represents the sign of `self`.
    ///
    /// - `1.0` if the number is positive, `+0.0` or `inf`
    /// - `-1.0` if the number is negative, `-0.0` or `-inf`
    /// - `NaN` if the number is `NaN`
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let f = FpDouble::from(3.5);
    ///
    /// assert_eq!(f.signum(), FpDouble::one());
    /// assert_eq!(FpDouble::neg_infinity().signum(), -FpDouble::one());
    ///
    /// assert!(FpDouble::nan().signum().is_nan());
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn signum(self) -> Self
    {
        if self.is_nan()
        {
            return self
        }
        Self::one().copysign(self)
    }

    /// Returns `true` if `self` has a positive sign, including `+0.0`, NaNs with
    /// positive sign bit and positive infinity. Note that IEEE 754 doesn't assign any
    /// meaning to the sign bit in case of a NaN, and as Rust doesn't guarantee that
    /// the bit pattern of NaNs are conserved over arithmetic operations, the result of
    /// `is_sign_positive` on a NaN might produce an unexpected result in some cases.
    /// See [explanation of NaN as a special value](f32) for more info.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let nan = FpDouble::nan();
    /// let neg_nan = -FpDouble::nan();
    ///
    /// let f = FpDouble::from(7.0);
    /// let g = FpDouble::from(-7.0);
    ///
    /// assert!(f.is_sign_positive());
    /// assert!(!g.is_sign_positive());
    /// assert!(nan.is_sign_positive());
    /// assert!(!neg_nan.is_sign_positive());
    /// ```
    #[must_use = "this returns the result of the operation, without modifying the original"]
    pub fn is_sign_positive(self) -> bool
    {
        self.sign_bit().is_zero()
    }

    /// Returns `true` if `self` has a negative sign, including `-0.0`, NaNs with
    /// negative sign bit and negative infinity. Note that IEEE 754 doesn't assign any
    /// meaning to the sign bit in case of a NaN, and as Rust doesn't guarantee that
    /// the bit pattern of NaNs are conserved over arithmetic operations, the result of
    /// `is_sign_negative` on a NaN might produce an unexpected result in some cases.
    /// See [explanation of NaN as a special value](f32) for more info.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let nan = FpDouble::nan();
    /// let neg_nan = -FpDouble::nan();
    ///
    /// let f = FpDouble::from(7.0);
    /// let g = FpDouble::from(-7.0);
    ///
    /// assert!(!f.is_sign_negative());
    /// assert!(g.is_sign_negative());
    /// assert!(!nan.is_sign_negative());
    /// assert!(neg_nan.is_sign_negative());
    /// ```
    #[must_use = "this returns the result of the operation, without modifying the original"]
    pub fn is_sign_negative(self) -> bool
    {
        !self.sign_bit().is_zero()
    }


    /// Returns the least number greater than `self`.
    ///
    /// Let `TINY` be the smallest representable positive value. Then,
    ///  - if `self.is_nan()`, this returns `self`;
    ///  - if `self` is [`NEG_INFINITY`], this returns [`MIN`];
    ///  - if `self` is `-TINY`, this returns -0.0;
    ///  - if `self` is -0.0 or +0.0, this returns `TINY`;
    ///  - if `self` is [`MAX`] or [`INFINITY`], this returns [`INFINITY`];
    ///  - otherwise the unique least value greater than `self` is returned.
    ///
    /// 
    /// The identity `x.next_up() == -(-x).next_down()` holds for all non-NaN `x`.
    /// When `x` is finite and the radix is 2, `x == x.next_up().next_down()` also holds.
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    /// 
    /// // epsilon is the difference between 1.0 and the next number up.
    /// assert_eq!(FpSingle::one().next_up(), FpSingle::one() + FpSingle::epsilon());
    /// // But not for most numbers.
    /// assert!(FpSingle::from(0.1).next_up() < FpSingle::from(0.1) + FpSingle::epsilon());
    /// assert_eq!(FpSingle::from(16777216.0).next_up(), FpSingle::from(16777218.0));
    /// ```
    ///
    /// [`NEG_INFINITY`]: Self::neg_infinity
    /// [`INFINITY`]: Self::infinity
    /// [`MIN`]: Self::min_value
    /// [`MAX`]: Self::max_value
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn next_up(self) -> Self
    {
        if self.is_nan()
        {
            return self
        }
        if self.is_infinite()
        {
            if self.is_sign_negative()
            {
                return Self::min_value()
            }
            return self
        }

        if self.is_zero()
        {
            return Self::from_bits(U::one())
        }

        let s = self.sign_bit();
        let mut e = self.exp_bits();
        let mut f = self.frac_bits();
        
        let is_subnormal = !(!e.is_zero() || !Self::IS_INT_IMPLICIT);
        
        if !is_subnormal //normal
        {
            f = f + (self.int_bits() << FRAC_SIZE);
        }

        let base = U::from(EXP_BASE).unwrap();

        if s.is_zero()
        {
            f = f + U::one();
            while f > U::one() << Self::MANTISSA_OP_SIZE
            {
                e = e + U::one();
                f = util::rounding_div(f, base);
            }
        }
        else
        {
            while e > U::zero() && f <= U::one() << Self::MANTISSA_OP_SIZE - Self::BASE_PADDING
            {
                e = e - U::one();
                f = f*base;
            }
            f = f - U::one();
        }
        
        if is_subnormal //subnormal
        {
            f = f << 1usize
        }

        if f.is_zero()
        {
            return if !s.is_zero() {-Self::zero()} else {Self::zero()}
        }

        if e.is_zero() && Self::IS_INT_IMPLICIT // subnormal
        {
            Fp::from_bits((s << Self::SIGN_POS) + util::rounding_div_2(f))
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

            Self::from_bits((s << Self::SIGN_POS) + (e << Self::EXP_POS) + f)
        }
    }
    
    /// Returns the greatest number less than `self`.
    ///
    /// Let `TINY` be the smallest representable positive value. Then,
    ///  - if `self.is_nan()`, this returns `self`;
    ///  - if `self` is [`INFINITY`], this returns [`MAX`];
    ///  - if `self` is `TINY`, this returns 0.0;
    ///  - if `self` is -0.0 or +0.0, this returns `-TINY`;
    ///  - if `self` is [`MIN`] or [`NEG_INFINITY`], this returns [`NEG_INFINITY`];
    ///  - otherwise the unique greatest value less than `self` is returned.
    ///
    /// The identity `x.next_down() == -(-x).next_up()` holds for all non-NaN `x`. When `x`
    /// is finite `x == x.next_down().next_up()` also holds.
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let x = FpSingle::one();
    /// // Clamp value into range [0, 1).
    /// let clamped = x.clamp(FpSingle::zero(), FpSingle::one().next_down());
    /// assert!(clamped < FpSingle::one());
    /// assert_eq!(clamped.next_up(), FpSingle::one());
    /// ```
    ///
    /// [`NEG_INFINITY`]: Self::neg_infinity
    /// [`INFINITY`]: Self::infinity
    /// [`MIN`]: Self::min_value
    /// [`MAX`]: Self::max_value
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn next_down(self) -> Self
    {
        if self.is_nan()
        {
            return self
        }
        if self.is_infinite()
        {
            if self.is_sign_positive()
            {
                return Self::max_value()
            }
            return self
        }

        if self.is_zero()
        {
            return Self::from_bits((U::one() << Self::SIGN_POS) + U::one())
        }

        let s = self.sign_bit();
        let mut e = self.exp_bits();
        let mut f = self.frac_bits();
        
        let is_subnormal = !(!e.is_zero() || !Self::IS_INT_IMPLICIT);
        
        if !is_subnormal //normal
        {
            f = f + (self.int_bits() << FRAC_SIZE);
        }

        let base = U::from(EXP_BASE).unwrap();

        if !s.is_zero()
        {
            f = f + U::one();
            while f > U::one() << Self::MANTISSA_OP_SIZE
            {
                e = e + U::one();
                f = util::rounding_div(f, base);
            }
        }
        else
        {
            while e > U::zero() && f <= U::one() << Self::MANTISSA_OP_SIZE - Self::BASE_PADDING
            {
                e = e - U::one();
                f = f*base;
            }
            f = f - U::one();
        }
        
        if is_subnormal //subnormal
        {
            f = f << 1usize
        }

        if f.is_zero()
        {
            return if !s.is_zero() {-Self::zero()} else {Self::zero()}
        }

        if e.is_zero() && Self::IS_INT_IMPLICIT // subnormal
        {
            Fp::from_bits((s << Self::SIGN_POS) + util::rounding_div_2(f))
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

            Fp::from_bits((s << Self::SIGN_POS) + (e << Self::EXP_POS) + f)
        }
    }
    
    /// Returns the maximum of the two numbers, propagating NaN.
    ///
    /// This returns NaN when *either* argument is NaN, as opposed to
    /// [`Fp::max`] which only returns NaN when *both* arguments are NaN.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    /// 
    /// let x = FpSingle::from(1.0);
    /// let y = FpSingle::from(2.0);
    ///
    /// assert_eq!(x.maximum(y), y);
    /// assert!(x.maximum(FpSingle::nan()).is_nan());
    /// ```
    ///
    /// If one of the arguments is NaN, then NaN is returned. Otherwise this returns the greater
    /// of the two numbers. For this operation, -0.0 is considered to be less than +0.0.
    /// Note that this follows the semantics specified in IEEE 754-2019.
    ///
    /// Also note that "propagation" of NaNs here doesn't necessarily mean that the bitpattern of a NaN
    /// operand is conserved; see [explanation of NaN as a special value](f32) for more info.
    #[must_use = "this returns the result of the comparison, without modifying either input"]
    #[inline]
    pub fn maximum(self, other: Self) -> Self
    {
        match (self.is_nan(), other.is_nan())
        {
            (true, true) => {
                let s1 = self.is_sign_negative();
                let s2 = other.is_sign_negative();
                if s1 == s2
                {
                    return self + other
                }
                return if !s1 {self} else {other}
            },
            (true, false) => return self + other,
            (false, true) => return self + other,
            (false, false) => ()
        }
        match self.total_cmp(other)
        {
            Ordering::Equal => if self.is_sign_positive() && other.is_sign_negative() { self } else { other },
            Ordering::Greater => self,
            Ordering::Less => other
        }
    }
    
    /// Returns the minimum of the two numbers, propagating NaN.
    ///
    /// This returns NaN when *either* argument is NaN, as opposed to
    /// [`Fp::min`] which only returns NaN when *both* arguments are NaN.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    /// 
    /// let x = FpSingle::from(1.0);
    /// let y = FpSingle::from(2.0);
    ///
    /// assert_eq!(x.minimum(y), x);
    /// assert!(x.minimum(FpSingle::nan()).is_nan());
    /// ```
    ///
    /// If one of the arguments is NaN, then NaN is returned. Otherwise this returns the lesser
    /// of the two numbers. For this operation, -0.0 is considered to be less than +0.0.
    /// Note that this follows the semantics specified in IEEE 754-2019.
    ///
    /// Also note that "propagation" of NaNs here doesn't necessarily mean that the bitpattern of a NaN
    /// operand is conserved; see [explanation of NaN as a special value](f32) for more info.
    #[must_use = "this returns the result of the comparison, without modifying either input"]
    #[inline]
    pub fn minimum(self, other: Self) -> Self
    {
        match (self.is_nan(), other.is_nan())
        {
            (true, true) => {
                let s1 = self.is_sign_negative();
                let s2 = other.is_sign_negative();
                if s1 == s2
                {
                    return self + other
                }
                return if s1 {self} else {other}
            },
            (true, false) => return self + other,
            (false, true) => return self + other,
            (false, false) => ()
        }
        match self.total_cmp(other)
        {
            Ordering::Equal => if self.is_sign_negative() && other.is_sign_positive() { self } else { other },
            Ordering::Greater => other,
            Ordering::Less => self
        }
    }
    
    /// Calculates the middle point of `self` and `rhs`.
    ///
    /// This returns NaN when *either* argument is NaN or if a combination of
    /// +inf and -inf is provided as arguments.
    ///
    /// # Examples
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    /// 
    /// assert_eq!(FpSingle::from(1.0).midpoint(FpSingle::from(4.0)), FpSingle::from(2.5));
    /// assert_eq!(FpSingle::from(-5.5).midpoint(FpSingle::from(8.0)), FpSingle::from(1.25));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn midpoint(self, other: Self) -> Self
    {
        let half = <Self as From<_>>::from(0.5);

        let lo = Self::min_positive_value()*<Self as From<_>>::from(2.0);
        let hi = Self::max_value()*half;

        let (a, b) = (self, other);
        let abs_a = a.abs();
        let abs_b = b.abs();

        if abs_a <= hi && abs_b <= hi
        {
            // Overflow is impossible
            (a + b)*half
        }
        else if abs_a < lo
        {
            // Not safe to halve a
            a + b*half
        }
        else if abs_b < lo
        {
            // Not safe to halve b
            a*half + b
        }
        else
        {
            // Not safe to halve a and b
            a*half + b*half
        }
    }

    /// Fused multiply-add. Computes `(self * a) + b`.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let m = FpDouble::from(10.0);
    /// let x = FpDouble::from(4.0);
    /// let b = FpDouble::from(60.0);
    ///
    /// // 100.0
    /// let abs_difference = (m.mul_add(x, b) - (m*x + b)).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn mul_add(self, a: Self, b: Self) -> Self
    {
        (self*a) + b
    }

    /// Take the reciprocal (inverse) of a number, `1/x`.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let x = FpDouble::from(2.0);
    /// let abs_difference = (x.recip() - (FpDouble::one()/x)).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn recip(self) -> Self
    {
        Self::one()/self
    }

    /// Raises a number to an integer power.
    ///
    /// Using this function is generally faster than using `powf`.
    /// It might have a different sequence of rounding operations than `powf`,
    /// so the results are not guaranteed to agree.
    /// 
    /// # Examples
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let x = FpDouble::from(2.0);
    /// let abs_difference = (x.powi(2) - x*x).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn powi<I: Int>(self, n: I) -> Self
    {
        util::powi(self, n)
    }
    
    /// Raise a number to an unsigned integer power.
    ///
    /// Using this function is generally faster than using `powf`.
    /// It might have a different sequence of rounding operations than `powf`,
    /// so the results are not guaranteed to agree.
    /// 
    /// # Examples
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let x = FpDouble::from(2.0);
    /// let abs_difference = (x.powu(2u32) - x*x).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn powu<I: UInt>(self, n: I) -> Self
    {
        util::powu(self, n)
    }

    /// Raises a number to a floating point power.
    ///
    /// # Examples
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let x = FpDouble::from(2.0);
    /// let abs_difference = (x.powf(FpDouble::from(2.0)) - x*x).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn powf(self, n: Self) -> Self
    {
        let xabs = self.abs();
        let nabs = n.abs();

        let exp_frac = U::from(util::count_digits_in_base(FRAC_SIZE + INT_SIZE, EXP_BASE)).unwrap();

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

    /// Returns the square root of a number.
    ///
    /// Returns NaN if `self` is a negative number other than `-0.0`.
    ///
    /// # Examples
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let positive = FpDouble::from(4.0);
    /// let negative = FpDouble::from(-4.0);
    ///
    /// let abs_difference = (positive.sqrt() - FpDouble::from(2.0)).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-10));
    /// assert!(negative.sqrt().is_nan());
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn sqrt(self) -> Self
    {   
        if self.is_nan()
        {
            return self
        }
        if self.is_zero()
        {
            return self.abs()
        }
        if self.is_sign_negative()
        {
            return Self::snan()
        }
        if self.is_infinite()
        {
            return self
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

        const NEWTON: usize = NEWTON_RT;
        let half = <Self as From<_>>::from(0.5);
        let mut y = y;
        for _ in 0..NEWTON
        {
            y = half*(y + self/y);
        }
        y
    }

    /// Returns `EXP_BASE^(self)`.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::{FpDouble, DecDouble};
    ///
    /// let f = FpDouble::from(2.0);
    /// let d = DecDouble::from(2.0);
    ///
    /// // 2^2 - 4 == 0
    /// let abs_difference_f = (f.exp_base() - FpDouble::from(4.0)).abs();
    ///
    /// // 10^2 - 100 == 0
    /// let abs_difference_d = (d.exp_base() - DecDouble::from(100.0)).abs();
    ///
    /// assert!(abs_difference_f < FpDouble::from(1e-10));
    /// assert!(abs_difference_d < DecDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
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
        
        let exp_frac = U::from(util::count_digits_in_base(FRAC_SIZE + INT_SIZE, EXP_BASE)).unwrap();

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
        let e = e.to_uint::<U>().unwrap();
        let i = if INT_SIZE > 0
        {
            U::one()
        }
        else
        {
            U::zero()
        };
        let y = z*Self::from_bits((e << Self::EXP_POS) + (i << Self::INT_POS));

        if neg
        {
            return y.recip()
        }
        y
    }

    fn exp_nonewton(self) -> Self
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

    /// Returns `e^(self)`, (the exponential function).
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let one = FpDouble::one();
    /// // e^1
    /// let e = one.exp();
    /// 
    /// // ln(e) - 1 == 0
    /// let abs_difference = (e.ln() - FpDouble::one()).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn exp(self) -> Self
    {
        let mut y = self.exp_nonewton();

        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_EXP;

            for _ in 0..NEWTON
            {
                y -= y*(y.ln_nonewton() - self)
            }
        }

        y
    }
    
    fn exp10_nonewton(self) -> Self
    {
        if EXP_BASE == 10
        {
            return self.exp_base()
        }
        (self*Self::LN_10()).exp_nonewton()
    }

    /// Returns `10^(self)`.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let f = FpDouble::from(2.0);
    ///
    /// // 10^2 - 100 == 0
    /// let abs_difference = (f.exp10() - FpDouble::from(100.0)).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-6));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn exp10(self) -> Self
    {
        if EXP_BASE == 10
        {
            return self.exp_base()
        }
        (self*Self::LN_10()).exp()
    }

    fn exp2_nonewton(self) -> Self
    {
        if EXP_BASE == 2
        {
            return self.exp_base()
        }
        (self*Self::LN_2()).exp_nonewton()
    }

    /// Returns `2^(self)`.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let f = FpDouble::from(2.0);
    ///
    /// // 2^2 - 4 == 0
    /// let abs_difference = (f.exp2() - FpDouble::from(4.0)).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn exp2(self) -> Self
    {
        if EXP_BASE == 2
        {
            return self.exp_base()
        }
        (self*Self::LN_2()).exp()
    }

    fn ln_nonewton(self) -> Self
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

    /// Returns the natural logarithm of the number.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let one = FpDouble::one();
    /// // e^1
    /// let e = one.exp();
    /// 
    /// // ln(e) - 1 == 0
    /// let abs_difference = (e.ln() - FpDouble::one()).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn ln(self) -> Self
    {
        let mut y = self.ln_nonewton();

        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_LN;

            for _ in 0..NEWTON
            {
                y -= Self::one() - self/y.exp_nonewton()
            }
        }

        y
    }

    /// Returns the logarithm of the number with respect to an arbitrary base.
    ///
    /// The result might not be correctly rounded owing to implementation details;
    /// `self.log2()` can produce more accurate results for base 2, and
    /// `self.log10()` can produce more accurate results for base 10.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let ten = FpDouble::from(10.0);
    /// let two = FpDouble::from(2.0);
    ///
    /// // log10(10) - 1 == 0
    /// let abs_difference_10 = (ten.log(FpDouble::from(10.0)) - FpDouble::one()).abs();
    ///
    /// // log2(2) - 1 == 0
    /// let abs_difference_2 = (two.log(FpDouble::from(2.0)) - FpDouble::one()).abs();
    ///
    /// assert!(abs_difference_10 < FpDouble::from(1e-10));
    /// assert!(abs_difference_2 < FpDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn log(self, base: Self) -> Self
    {
        self.ln()/base.ln()
    }

    /// Returns the logarithm base `EXP_BASE` of the number.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::{FpDouble, DecDouble};
    ///
    /// let two = FpDouble::from(2.0);
    /// let ten = DecDouble::from(10.0);
    ///
    /// // log2(2) - 1 == 0
    /// let abs_difference_2 = (two.log_base() - FpDouble::one()).abs();
    /// 
    /// // log10(10) - 1 == 0
    /// let abs_difference_10 = (ten.log_base() - DecDouble::one()).abs();
    ///
    /// assert!(abs_difference_2 < FpDouble::from(1e-10));
    /// assert!(abs_difference_10 < DecDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
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
    
    fn log2_nonewton(self) -> Self
    {
        if EXP_BASE == 2
        {
            return self.log_base()
        }
        self.ln_nonewton()/Self::LN_2()
    }

    /// Returns the base 2 logarithm of the number.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let two = FpDouble::from(2.0);
    ///
    /// // log2(2) - 1 == 0
    /// let abs_difference = (two.log2() - FpDouble::one()).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn log2(self) -> Self
    {
        let mut y = self.log2_nonewton();

        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_LN;

            for _ in 0..NEWTON
            {
                y -= (Self::one() - self/y.exp2_nonewton())/Self::LN_2()
            }
        }

        y
    }
    
    fn log10_nonewton(self) -> Self
    {
        if EXP_BASE == 10
        {
            return self.log_base()
        }
        self.ln_nonewton()/Self::LN_10()
    }

    /// Returns the base 10 logarithm of the number.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let ten = FpDouble::from(10.0);
    ///
    /// // log10(10) - 1 == 0
    /// let abs_difference = (ten.log10() - FpDouble::one()).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-3));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn log10(self) -> Self
    {
        let mut y = self.log10_nonewton();

        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_LN;

            for _ in 0..NEWTON
            {
                y -= (Self::one() - self/y.exp10_nonewton())/Self::LN_10()
            }
        }

        y
    }

    /// Converts radians to degrees.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use num::traits::FloatConst;
    ///
    /// let angle = FpDouble::PI();
    ///
    /// let abs_difference = (angle.to_degrees() - FpDouble::from(180.0)).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn to_degrees(self) -> Self
    {
        self/(Self::FRAC_PI_2()/<Self as From<_>>::from(90.0))
    }
    
    /// Converts degrees to radians.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use num::traits::FloatConst;
    ///
    /// let angle = FpDouble::from(180.0);
    ///
    /// let abs_difference = (angle.to_radians() - FpDouble::PI()).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn to_radians(self) -> Self
    {
        self*(Self::FRAC_PI_2()/<Self as From<_>>::from(90.0))
    }

    /// Returns the maximum of the two numbers.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let x = FpDouble::from(1.0);
    /// let y = FpDouble::from(2.0);
    ///
    /// assert_eq!(x.max(y), y);
    /// ```
    #[must_use = "this returns the result of the comparison, without modifying either input"]
    pub fn max(self, other: Self) -> Self
    {
        match (self.is_nan(), other.is_nan())
        {
            (true, true) => {
                let s1 = self.is_sign_negative();
                let s2 = other.is_sign_negative();
                if s1 == s2
                {
                    return self + other
                }
                return if !s1 {self} else {other}
            },
            (true, false) => return other,
            (false, true) => return self,
            (false, false) => ()
        }
        match self.total_cmp(other)
        {
            Ordering::Equal => if self.is_sign_positive() && other.is_sign_negative() { self } else { other },
            Ordering::Greater => self,
            Ordering::Less => other
        }
    }

    /// Returns the minimum of the two numbers.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let x = FpDouble::from(1.0);
    /// let y = FpDouble::from(2.0);
    ///
    /// assert_eq!(x.min(y), x);
    /// ```
    #[must_use = "this returns the result of the comparison, without modifying either input"]
    pub fn min(self, other: Self) -> Self
    {
        match (self.is_nan(), other.is_nan())
        {
            (true, true) => {
                let s1 = self.is_sign_negative();
                let s2 = other.is_sign_negative();
                if s1 == s2
                {
                    return self + other
                }
                return if s1 {self} else {other}
            },
            (true, false) => return other,
            (false, true) => return self,
            (false, false) => ()
        }
        match self.total_cmp(other)
        {
            Ordering::Equal => if self.is_sign_negative() && other.is_sign_positive() { self } else { other },
            Ordering::Greater => other,
            Ordering::Less => self
        }
    }

    /// The positive difference of two numbers.
    ///
    /// * If `self <= other`: `0:0`
    /// * Else: `self - other`
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let x = FpDouble::from(3.0);
    /// let y = FpDouble::from(-3.0);
    ///
    /// let abs_difference_x = (x.abs_sub(FpDouble::one()) - FpDouble::from(2.0)).abs();
    /// let abs_difference_y = (y.abs_sub(FpDouble::one()) - FpDouble::zero()).abs();
    ///
    /// assert!(abs_difference_x < FpDouble::from(1e-10));
    /// assert!(abs_difference_y < FpDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[deprecated(
        note = "you probably meant `(self - other).abs()`: \
                this operation is `(self - other).max(0.0)` \
                except that `abs_sub` also propagates NaNs (also \
                known as `fdimf` in C). If you truly need the positive \
                difference, consider using that expression or the C function \
                `fdimf`, depending on how you wish to handle NaN (please consider \
                filing an issue describing your use-case too)."
    )]
    pub fn abs_sub(self, other: Self) -> Self
    {
        (self - other).max(Self::zero())
    }

    /// Take the cubic root of a number.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let x = FpDouble::from(8.0);
    ///
    /// // x^(1/3) - 2 == 0
    /// let abs_difference = (x.cbrt() - FpDouble::from(2.0)).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-9));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn cbrt(self) -> Self
    {
        if self.is_nan()
        {
            return self
        }
        if self.is_infinite()
        {
            return self
        }
        if self.is_zero()
        {
            return self
        }

        let y = {
            let xabs_log = self.abs().log_base();

            let n_xabs_log = xabs_log/<Self as From<_>>::from(3.0);
            
            n_xabs_log.exp_base().copysign(self)
        };

        const NEWTON: usize = NEWTON_RT;
        let third = <Self as From<_>>::from(3.0).recip();
        let two = <Self as From<_>>::from(2.0);
        let mut y = y;
        for _ in 0..NEWTON
        {
            y = third*(self/(y*y) + two*y);
        }
        y
    }

    /// Compute the distance between the origin and a point (`x`, `y`) on the
    /// Euclidean plane. Equivalently, compute the length of the hypotenuse of a
    /// right-angle triangle with other sides having length `x.abs()` and
    /// `y.abs()`.
    ///
    /// # Examples
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let x = FpDouble::from(2.0);
    /// let y = FpDouble::from(3.0);
    ///
    /// // sqrt(x^2 + y^2)
    /// let abs_difference = (x.hypot(y) - (x.powi(2) + y.powi(2)).sqrt()).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn hypot(self, other: Self) -> Self
    {
        (self*self + other*other).sqrt()
    }

    /// Computes the sine of a number (in radians).
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use num::traits::FloatConst;
    ///
    /// let x = FpDouble::FRAC_PI_2();
    ///
    /// let abs_difference = (x.sin() - FpDouble::one()).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
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

        static mut P: Option<[f64; N]> = None;
        let p = unsafe {
            if P.is_none()
            {
                P = Some({
                    let t: [[f64; N]; N] = ArrayOps::fill(
                        |n| <[_; N]>::chebyshev_polynomial(1, n).unwrap()
                    );
                    let p: [f64; N] = t.zip(C)
                        .map2(|(t, c)| t.map2(|tn| c*tn))
                        .reduce(|a, b| a.zip(b).map2(|(a, b)| a + b))
                        .unwrap_or_default();
                    p
                })
            }
            P.unwrap()
        }.map(|p| <Self as From<_>>::from(p));

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

    /// Computes the cosine of a number (in radians).
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use num::traits::FloatConst;
    ///
    /// let x = FpDouble::TAU();
    ///
    /// let abs_difference = (x.cos() - FpDouble::one()).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-9));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
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

        static mut P: Option<[f64; N]> = None;
        let p = unsafe {
            if P.is_none()
            {
                P = Some({
                    let t: [[f64; N]; N] = ArrayOps::fill(
                        |n| <[_; N]>::chebyshev_polynomial(1, n).unwrap()
                    );
                    let p: [f64; N] = t.zip(C)
                        .map2(|(t, c)| t.map2(|tn| c*tn))
                        .reduce(|a, b| a.zip(b).map2(|(a, b)| a + b))
                        .unwrap_or_default();
                    p
                })
            }
            P.unwrap()
        }.map(|p| <Self as From<_>>::from(p));

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

    /// Computes the tangent of a number (in radians).
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use num::traits::FloatConst;
    ///
    /// let x = FpDouble::FRAC_PI_4();
    /// 
    /// let abs_difference = (x.tan() - FpDouble::one()).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-3));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn tan(self) -> Self
    {
        let (sin, cos) = self.sin_cos();

        sin/cos
    }

    /// Computes the arcsine of a number. Return value is in radians in
    /// the range [-pi/2, pi/2] or NaN if the number is outside the range
    /// [-1, 1].
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use num::traits::FloatConst;
    ///
    /// let f = FpDouble::FRAC_PI_2();
    ///
    /// // asin(sin(pi/2))
    /// let abs_difference = (f.sin().asin() - f).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-5));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
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
        let xabs = self.abs();
        if xabs > Self::one()
        {
            return (self - self)/(self - self)
        }

        let mut y = (self/(Self::one() - self*self).sqrt()).atan();

        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_TRIG;

            for _ in 0..NEWTON
            {
                let (sin, cos) = y.sin_cos();
                y -= (sin - self)/cos
            }
        }

        y
    }

    /// Computes the arccosine of a number. Return value is in radians in
    /// the range [0, pi] or NaN if the number is outside the range
    /// [-1, 1].
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use num::traits::FloatConst;
    ///
    /// let f = FpDouble::FRAC_PI_4();
    ///
    /// // acos(cos(pi/4))
    /// let abs_difference = (f.cos().acos() - f).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-9));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
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
        let xabs = self.abs();
        if xabs > Self::one()
        {
            return (self - self)/(self - self)
        }
        let mut y = ((Self::one() - self*self).sqrt()/self).atan();

        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_TRIG;

            for _ in 0..NEWTON
            {
                let (sin, cos) = y.sin_cos();
                y += (cos - self)/sin
            }
        }

        y
    }

    /// Computes the arctangent of a number. Return value is in radians in the
    /// range [-pi/2, pi/2];
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let f = FpDouble::one();
    ///
    /// // atan(tan(1))
    /// let abs_difference = (f.tan().atan() - f).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-4));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
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
        let mut y = if self.abs() < Self::one()
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
        };

        y = y % Self::PI();
        while y > Self::FRAC_PI_2()
        {
            y -= Self::PI()
        }
        while y < -Self::FRAC_PI_2()
        {
            y += Self::PI()
        }
        
        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_TRIG;

            for _ in 0..NEWTON
            {
                let (sin, cos) = y.sin_cos();
                y -= (sin - cos*self)*cos
            }
        }
        
        if y.abs() > Self::FRAC_PI_2()
        {
            return Self::FRAC_PI_2().copysign(y)
        }

        y
    }

    /// Computes the four quadrant arctangent of `self` (`y`) and `other` (`x`).
    ///
    /// * `x = 0`, `y = 0`: `0`
    /// * `x >= 0`: `arctan(y/x)` -> `[-pi/2, pi/2]`
    /// * `y >= 0`: `arctan(y/x) + pi` -> `(pi/2, pi]`
    /// * `y < 0`: `arctan(y/x) - pi` -> `(-pi, -pi/2)`
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use num::traits::FloatConst;
    ///
    /// // All angles from horizontal right (+x)
    /// // 45 deg counter-clockwise
    /// let x1 = FpDouble::from(3.0);
    /// let y1 = FpDouble::from(-3.0);
    ///
    /// // 135 deg clockwise
    /// let x2 = FpDouble::from(-3.0);
    /// let y2 = FpDouble::from(3.0);
    ///
    /// let abs_difference_1 = (y1.atan2(x1) - (-FpDouble::FRAC_PI_4())).abs();
    /// let abs_difference_2 = (y2.atan2(x2) - (FpDouble::PI() - FpDouble::FRAC_PI_4())).abs();
    ///
    /// assert!(abs_difference_1 < FpDouble::from(1e-10));
    /// assert!(abs_difference_2 < FpDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn atan2(self, other: Self) -> Self
    {
        if other.is_zero()
        {
            return if self.is_zero()
            {
                Self::zero()
            }
            else if self >= Self::zero()
            {
                Self::FRAC_PI_2()
            }
            else
            {
                -Self::FRAC_PI_2()
            }
        }
        let atan = (self/other).atan();
        if other >= Self::zero()
        {
            return atan
        }
        atan + if self >= Self::zero()
        {
            Self::PI()
        }
        else
        {
            -Self::PI()
        }
    }

    /// Simultaneously computes the sine and cosine of the number, `x`. Returns `(sin(x), cos(x))`.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use num::traits::FloatConst;
    ///
    /// let x = FpDouble::FRAC_PI_4();
    /// let f = x.sin_cos();
    ///
    /// let abs_difference_0 = (f.0 - x.sin()).abs();
    /// let abs_difference_1 = (f.1 - x.cos()).abs();
    ///
    /// assert!(abs_difference_0 < FpDouble::from(1e-10));
    /// assert!(abs_difference_1 < FpDouble::from(1e-10));
    /// ```
    #[must_use]
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

        static mut P: Option<[f64; N]> = None;
        let p = unsafe {
            if P.is_none()
            {
                P = Some({
                    let t: [[f64; N]; N] = ArrayOps::fill(
                        |n| <[_; N]>::chebyshev_polynomial(1, n).unwrap()
                    );
                    let p: [f64; N] = t.zip(C)
                        .map2(|(t, c)| t.map2(|tn| c*tn))
                        .reduce(|a, b| a.zip(b).map2(|(a, b)| a + b))
                        .unwrap_or_default();
                    p
                })
            }
            P.unwrap()
        }.map(|p| <Self as From<_>>::from(p));

        let two = <Self as From<_>>::from(2.0);

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
            let w = if w > Self::one() {two - w} else if w < -Self::one() {-two - w} else {w};

            let z = two*w*w - Self::one();

            p.polynomial(z)*w
        };

        (sin, cos)
    }

    /// Returns `e^(self) - 1`.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let x = FpDouble::from(7.0);
    ///
    /// // e^(ln(7)) - 1
    /// let abs_difference = (x.ln().exp_m1() - FpDouble::from(6.0)).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-8));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn exp_m1(self) -> Self
    {
        let mut y = self.exp_nonewton() - Self::one();

        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_EXP;

            for _ in 0..NEWTON
            {
                let yp1 = y + Self::one();
                y -= yp1*(yp1.ln() - self)
            }
        }

        y
    }

    /// Returns `ln(1+n)` (natural logarithm) more accurately than if
    /// the operations were performed separately.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use num::traits::FloatConst;
    ///
    /// let x = FpDouble::E() - FpDouble::one();
    ///
    /// // ln(1 + (e - 1)) == ln(e) == 1
    /// let abs_difference = (x.ln_1p() - FpDouble::one()).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-9));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn ln_1p(self) -> Self
    {
        let xp1 = self + Self::one();
        let mut y = xp1.ln_nonewton();

        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_LN;

            for _ in 0..NEWTON
            {
                y -= Self::one() - xp1/y.exp_nonewton()
            }
        }

        y
    }

    /// Hyperbolic sine function.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use num::traits::FloatConst;
    ///
    /// let e = FpDouble::E();
    /// let x = FpDouble::one();
    ///
    /// let f = x.sinh();
    /// // Solving sinh() at 1 gives `(e^2-1)/(2e)`
    /// let g = (e*e - FpDouble::one())/(FpDouble::from(2.0)*e);
    /// let abs_difference = (f - g).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn sinh(self) -> Self
    {
        if !self.is_finite()
        {
            return self
        }

        let emx = (-self.abs()).exp();
    
        ((Self::one() - emx*emx)/emx*<Self as From<_>>::from(0.5)).copysign(self)
    }

    /// Hyperbolic cosine function.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use num::traits::FloatConst;
    ///
    /// let e = FpDouble::E();
    /// let x = FpDouble::one();
    /// 
    /// let f = x.cosh();
    /// // Solving cosh() at 1 gives this result
    /// let g = (e*e + FpDouble::one())/(FpDouble::from(2.0)*e);
    /// let abs_difference = (f - g).abs();
    ///
    /// // Same result
    /// assert!(abs_difference < FpDouble::from(1.0e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn cosh(self) -> Self
    {
        if !self.is_finite()
        {
            return self.abs()
        }

        let emx = (-self.abs()).exp();
    
        (Self::one() + emx*emx)/emx*<Self as From<_>>::from(0.5)
    }
    
    /// Hyperbolic tangent function.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use num::traits::FloatConst;
    ///
    /// let e = FpDouble::E();
    /// let x = FpDouble::one();
    ///
    /// let f = x.tanh();
    /// // Solving tanh() at 1 gives `(1 - e^(-2))/(1 + e^(-2))`
    /// let g = (FpDouble::one() - e.powi(-2))/(FpDouble::one() + e.powi(-2));
    /// let abs_difference = (f - g).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1.0e-3));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn tanh(self) -> Self
    {
        if self.is_nan()
        {
            return self
        }
        if self.is_infinite()
        {
            return Self::one().copysign(self)
        }

        let ex = (-self.abs()).exp();
        let ex2 = ex*ex;
        let ex2p1 = Self::one() + ex2;
        let ex2m1 = Self::one() - ex2;
        
        (ex2m1/ex2p1).copysign(self)
    }

    /// Inverse hyperbolic sine function.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let x = FpDouble::one();
    /// let f = x.sinh().asinh();
    ///
    /// let abs_difference = (f - x).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1.0e-3));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn asinh(self) -> Self
    {
        if !self.is_finite()
        {
            return self
        }
        (self + (self*self + Self::one()).sqrt()).ln()
    }

    /// Inverse hyperbolic cosine function.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let x = FpDouble::one();
    /// let f = x.cosh().acosh();
    ///
    /// let abs_difference = (f - x).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1.0e-4));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn acosh(self) -> Self
    {
        if self.is_nan()
        {
            return self
        }
        if self < Self::one()
        {
            return (self - self)/(self - self)
        }
        if self.is_infinite()
        {
            return self
        }
        (self + (self*self - Self::one()).sqrt()).ln()
    }

    /// Inverse hyperbolic tangent function.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use num::traits::FloatConst;
    ///
    /// let e = FpDouble::E();
    /// let f = e.tanh().atanh();
    ///
    /// let abs_difference = (f - e).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1.0e-2));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
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
        <Self as From<_>>::from(0.5)*((Self::one() + self.abs())/(Self::one() - self.abs())).ln().copysign(self)
    }
    
    fn ln_gamma_lanczos(self) -> Self
    {
        let xm1 = self - Self::one();

        const LANCZOS_CHEB_7: [f64; 9] = [
            0.99999999999980993227684700473478,
            676.520368121885098567009190444019,
            -1259.13921672240287047156078755283,
            771.3234287776530788486528258894,
            -176.61502916214059906584551354,
            12.507343278686904814458936853,
            -0.13857109526572011689554707,
            9.984369578019570859563e-6,
            1.50563273514931155834e-7
        ];

        const LOGROOT2PI: f64 = 0.9189385332046727417803297364056176;

        let mut sum = <Self as From<_>>::from(LANCZOS_CHEB_7[0]);
        for k in 1..LANCZOS_CHEB_7.len()
        {
            sum += <Self as From<_>>::from(LANCZOS_CHEB_7[k])/(xm1 + Self::from_uint(k))
        }

        let term1 = (xm1 + <Self as From<_>>::from(0.5))
            *((xm1 + <Self as From<_>>::from(7.5))/Self::E()).ln();
        let term2 = <Self as From<_>>::from(LOGROOT2PI) + sum.ln();

        term1 + (term2 - <Self as From<_>>::from(7u8))
    }

    /// Natural logarithm of the absolute value of the gamma function
    ///
    /// The integer part of the tuple indicates the sign of the gamma function.
    ///
    /// # Examples
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// 
    /// let x = FpDouble::from(2.0);
    ///
    /// let abs_difference = (x.ln_gamma().0 - FpDouble::zero()).abs();
    ///
    /// assert!(abs_difference <= FpDouble::from(1e-2));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn ln_gamma(self) -> (Self, i32)
    {
        if self > <Self as From<_>>::from(0.5)
        {
            return (self.ln_gamma_lanczos(), 1)
        }

        let sin_fact = (Self::PI()*self).sin();

        if sin_fact.is_zero()
        {
            return (Self::snan(), 0)
        }

        (
            Self::PI().ln() - sin_fact.abs().ln() + (Self::one() - self).ln_gamma_lanczos(),
            if sin_fact.is_sign_negative() {-1} else {1}
        )
    }

    /// Gamma function.
    ///
    /// # Examples
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// 
    /// let x = FpDouble::from(5.0f32);
    ///
    /// let abs_difference = (x.gamma() - FpDouble::from(24.0)).abs();
    ///
    /// assert!(abs_difference <= FpDouble::from(1e-1));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn gamma(self) -> Self
    {
        let (lgamma, sign) = self.ln_gamma();
        if sign == 0
        {
            return Self::snan()
        }

        let gamma = lgamma.exp();

        if sign < 0
        {
            return -gamma
        }
        gamma
    }

    /// Returns a number composed of the magnitude of `self` and the sign of
    /// `sign`.
    ///
    /// Equal to `self` if the sign of `self` and `sign` are the same, otherwise
    /// equal to `-self`. If `self` is a NaN, then a NaN with the sign bit of
    /// `sign` is returned. Note, however, that conserving the sign bit on NaN
    /// across arithmetical operations is not generally guaranteed.
    ///
    /// # Examples
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let f = FpDouble::from(3.5);
    /// let s = FpDouble::from(0.42);
    ///
    /// assert_eq!(f.copysign(s), f);
    /// assert_eq!(f.copysign(-s), -f);
    /// assert_eq!((-f).copysign(s), f);
    /// assert_eq!((-f).copysign(-s), -f);
    ///
    /// assert!(FpDouble::nan().copysign(FpDouble::one()).is_nan());
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn copysign(self, sign: Self) -> Self
    {
        let mask = (U::max_value() >> util::bitsize_of::<U>() - Self::SIGN_SIZE) << Self::SIGN_POS;
        
        Self::from_bits((self.to_bits() & (!mask)) | (sign.to_bits() & mask))
    }
    
    /// Returns the additive identity element of `Self`, `0`.
    #[must_use]
    pub fn zero() -> Self
    {
        Self::from_bits(U::zero())
    }
    
    /// Sets `self` to the additive identity element of `Self`, `0`.
    fn set_zero(&mut self)
    {
        *self = Self::zero()
    }

    /// Returns `true` if `self` is equal to the additive identity.
    #[must_use = "this returns the result of the operation, without modifying the original"]
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
    
    /// Returns the multiplicative identity element of `Self`, `1`.
    #[must_use]
    pub fn one() -> Self
    {
        let bias = Self::exp_bias();
        if !Self::IS_INT_IMPLICIT
        {
            return Self::from_bits((bias << Self::EXP_POS) + (U::one() << Self::INT_POS))
        }
        Self::from_bits(bias << Self::EXP_POS)
    }
    
    /// Sets `self` to the multiplicative identity element of `Self`, `1`.
    fn set_one(&mut self)
    {
        *self = Self::one()
    }

    /// Returns `true` if `self` is equal to the multiplicative identity.
    #[must_use = "this returns the result of the operation, without modifying the original"]
    pub fn is_one(self) -> bool
    {
        self == Self::one()
    }
    
    /// Calculates Euclidean division, the matching method for `rem_euclid`.
    ///
    /// This computes the integer `n` such that
    /// `self = n * rhs + self.rem_euclid(rhs)`.
    /// In other words, the result is `self / rhs` rounded to the integer `n`
    /// such that `self >= n * rhs`.
    ///
    /// # Examples
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// 
    /// let a = FpDouble::from(7.0);
    /// let b = FpDouble::from(4.0);
    /// assert_eq!(a.div_euclid(b), FpDouble::from(1.0)); // 7.0 > 4.0 * 1.0
    /// assert_eq!((-a).div_euclid(b), FpDouble::from(-2.0)); // -7.0 >= 4.0 * -2.0
    /// assert_eq!(a.div_euclid(-b), FpDouble::from(-1.0)); // 7.0 >= -4.0 * -1.0
    /// assert_eq!((-a).div_euclid(-b), FpDouble::from(2.0)); // -7.0 >= -4.0 * 2.0
    /// ```
    #[doc(alias = "modulo", alias = "mod")]
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn div_euclid(self, rhs: Self) -> Self
    {
        let q = (self / rhs).trunc();
        if self % rhs < Self::zero()
        {
            return if rhs > Self::zero() {q - Self::one()} else {q + Self::one()};
        }
        q
    }

    /// Calculates the least nonnegative remainder of `self (mod rhs)`.
    ///
    /// In particular, the return value `r` satisfies `0.0 <= r < rhs.abs()` in
    /// most cases. However, due to a floating point round-off error it can
    /// result in `r == rhs.abs()`, violating the mathematical definition, if
    /// `self` is much smaller than `rhs.abs()` in magnitude and `self < 0.0`.
    /// This result is not an element of the function's codomain, but it is the
    /// closest floating point number in the real numbers and thus fulfills the
    /// property `self == self.div_euclid(rhs) * rhs + self.rem_euclid(rhs)`
    /// approximately.
    ///
    /// # Examples
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// 
    /// let a = FpDouble::from(7.0);
    /// let b = FpDouble::from(4.0);
    /// assert_eq!(a.rem_euclid(b), FpDouble::from(3.0));
    /// assert_eq!((-a).rem_euclid(b), FpDouble::from(1.0));
    /// assert_eq!(a.rem_euclid(-b), FpDouble::from(3.0));
    /// assert_eq!((-a).rem_euclid(-b), FpDouble::from(1.0));
    /// // limitation due to round-off error
    /// assert!((-FpDouble::epsilon()).rem_euclid(FpDouble::from(3.0)) != FpDouble::from(0.0));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn rem_euclid(self, rhs: Self) -> Self
    {
        let r = self % rhs;
        if r < Self::zero() {r + rhs.abs()} else {r}
    }
    
    /// Return the ordering between `self` and `other`.
    ///
    /// Unlike the standard partial comparison between floating point numbers,
    /// this comparison always produces an ordering in accordance to
    /// the `totalOrder` predicate as defined in the IEEE 754 (2008 revision)
    /// floating point standard. The values are ordered in the following sequence:
    ///
    /// - negative quiet NaN
    /// - negative signaling NaN
    /// - negative infinity
    /// - negative numbers
    /// - negative subnormal numbers
    /// - negative zero
    /// - positive zero
    /// - positive subnormal numbers
    /// - positive numbers
    /// - positive infinity
    /// - positive signaling NaN
    /// - positive quiet NaN.
    ///
    /// The ordering established by this function does not always agree with the
    /// [`PartialOrd`] and [`PartialEq`] implementations. For example,
    /// they consider negative and positive zero equal, while `total_cmp`
    /// doesn't.
    ///
    /// The interpretation of the signaling NaN bit follows the definition in
    /// the IEEE 754 standard, which may not match the interpretation by some of
    /// the older, non-conformant (e.g. MIPS) hardware implementations.
    ///
    /// # Examples
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::{FpSingle, FpDouble};
    /// use std::cmp::Ordering;
    ///
    /// assert_eq!(FpDouble::nan().total_cmp(FpDouble::nan()), Ordering::Equal);
    /// assert_eq!(FpSingle::nan().total_cmp(FpSingle::nan()), Ordering::Equal);
    ///
    /// assert_eq!((-FpDouble::nan()).total_cmp(FpDouble::nan()), Ordering::Less);
    /// assert_eq!(FpDouble::infinity().total_cmp(FpDouble::nan()), Ordering::Less);
    /// assert_eq!((-FpDouble::zero()).total_cmp(FpDouble::zero()), Ordering::Less);
    /// ```
    #[must_use = "this returns the result of the comparison, without modifying either input"]
    pub fn total_cmp(self, other: Self) -> std::cmp::Ordering
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
            
            while e0 > e1 && f0 < U::one() << Self::MANTISSA_OP_SIZE - Self::BASE_PADDING
            {
                e0 = e0 - U::one();
                f0 = f0*base;
            }            
            while e1 > e0 && f1 < U::one() << Self::MANTISSA_OP_SIZE - Self::BASE_PADDING
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

    /// Returns `self*EXP_BASE`.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::{FpDouble, DecDouble};
    ///
    /// let f = FpDouble::from(2.0);
    /// let d = DecDouble::from(2.0);
    ///
    /// // 2*2 - 4 == 0
    /// let abs_difference_f = (f.mul_base() - FpDouble::from(4.0)).abs();
    ///
    /// // 2*10 - 20 == 0
    /// let abs_difference_d = (d.mul_base() - DecDouble::from(20.0)).abs();
    ///
    /// assert!(abs_difference_f < FpDouble::from(1e-10));
    /// assert!(abs_difference_d < DecDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
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
        while e > U::zero() && f < U::one() << Self::MANTISSA_OP_SIZE - Self::BASE_PADDING
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

    /// Returns `self/EXP_BASE`.
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::{FpDouble, DecDouble};
    ///
    /// let f = FpDouble::from(2.0);
    /// let d = DecDouble::from(2.0);
    ///
    /// // 2/2 - 1 == 0
    /// let abs_difference_f = (f.div_base() - FpDouble::from(1.0)).abs();
    ///
    /// // 2/10 - 0.2 == 0
    /// let abs_difference_d = (d.div_base() - DecDouble::from(0.2)).abs();
    ///
    /// assert!(abs_difference_f < FpDouble::from(1e-10));
    /// assert!(abs_difference_d < DecDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
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
              
        while e > U::zero() && f < U::one() << Self::MANTISSA_OP_SIZE - Self::BASE_PADDING
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

    /// Restrict a value to a certain interval unless it is NaN.
    ///
    /// Returns `max` if `self` is greater than `max`, and `min` if `self` is
    /// less than `min`. Otherwise this returns `self`.
    ///
    /// Note that this function returns NaN if the initial value was NaN as
    /// well.
    ///
    /// # Panics
    ///
    /// Panics if `min > max`, `min` is NaN, or `max` is NaN.
    ///
    /// # Examples
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let min = FpSingle::from(-2.0);
    /// let max = FpSingle::one();
    ///
    /// assert_eq!(FpSingle::from(-3.0).clamp(min, max), FpSingle::from(-2.0));
    /// assert_eq!(FpSingle::zero().clamp(min, max), FpSingle::zero());
    /// assert_eq!(FpSingle::from(2.0).clamp(min, max), FpSingle::one());
    /// assert!(FpSingle::nan().clamp(min, max).is_nan());
    /// ```
    #[must_use = "this returns the result of the comparison, without modifying either input"]
    pub fn clamp(mut self, min: Self, max: Self) -> Self
    {
        assert!(min <= max, "min > max, or either was NaN. min = {min:?}, max = {max:?}");
        if self < min
        {
            self = min;
        }
        if self > max
        {
            self = max;
        }
        self
    }
}

impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:,
    [(); EXP_BASE - 2]:,
    U: ConstZero
{
    /// The additive identity element of `Self`, `0`.
    pub const ZERO: Self = Self::from_bits(U::ZERO);
}

#[cfg(test)]
mod test
{
    use crate::ieee754::{FpDouble, FpHalf};

    #[test]
    fn test_gamma()
    {
        let x = FpDouble::from(5.0f32);
        
        println!("{}", x.gamma())
    }

    #[test]
    fn test_next_up_down()
    {
        type F = FpHalf;

        let mut x = F::min_value();

        while x < F::max_value()
        {
            let y = x.next_up();
            assert_ne!(x, y);
            assert_eq!(-(-x).next_down(), y);
            assert_eq!(x, y.next_down());
            x = y;
        }

        while x > F::min_value()
        {
            let y = x.next_down();
            assert_ne!(x, y);
            assert_eq!(-(-x).next_up(), y);
            assert_eq!(x, y.next_up());
            x = y;
        }
    }
}