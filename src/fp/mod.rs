use core::cmp::Ordering;
use core::num::FpCategory;
use core::ops::Neg;

use num_traits::{Zero, NumCast, ConstZero, FloatConst, FromBytes, ToBytes};

use crate::{util, AnyInt, Int, UInt};

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

const NO_NEWTON: bool = false;
const NEWTON_EXP: usize = if NO_NEWTON {0} else {3};
const NEWTON_LN: usize = if NO_NEWTON {0} else {2};
const NEWTON_RT: usize = if NO_NEWTON {0} else {4};
const NEWTON_TRIG: usize = if NO_NEWTON {0} else {3};

macro_rules! as_lossless {
    ($value:expr, $fn_as_lossless:expr, $fn:block) => {
        {
            #[cfg(any(test, feature = "use_std_float"))]
            let _as_lossless = crate::Fp::as_lossless(
                $value,
                $fn_as_lossless,
                $fn_as_lossless,
                $fn_as_lossless,
                $fn_as_lossless
            );
            #[cfg(all(not(test), feature = "use_std_float"))]
            if let Some([as_lossless]) = _as_lossless
            {
                return as_lossless
            }
            let y = (|| $fn)();
            #[cfg(test)]
            if let Some([as_lossless]) = _as_lossless
            {
                if !y.approx_eq(as_lossless)
                {
                    debug_assert_eq!(y, as_lossless, "Error is too big!")
                }
            }
            y
        }
    };
}
use as_lossless as as_lossless;

/// A custom floating point type, where the bit size of the exponent and mantissa can be set separately.
/// 
/// `U` is the underlying unsigned integer type which is used to represent the number.
///
/// `SIGN_BIT` is wether or not the number has a sign bit.
/// 
/// `EXP_SIZE` is the size of the exponent in bits.
/// 
/// `INT_SIZE` is the size of the integer part of the mantissa in bits. If zero, then the integer bit is implicit.
/// 
/// `FRAC_SIZE` is the size of the fractional part of the mantissa in bits.
/// 
/// `EXP_BASE` is the base of the exponent.
/// 
/// The total bit size of `U` must be greater or equal to `SIGN_BIT` + `EXP_SIZE` + `INT_SIZE` + `FRAC_SIZE` to contain the entire number.
///
/// The bit layout is as follows:
/// ```txt
/// No data: | Sign:      | Exponent:  | Integer:   | Fractional: |
/// <  ..  > | <SIGN_BIT> | <EXP_SIZE> | <INT_SIZE> | <FRAC_SIZE> |
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
/// 
/// The number then automatically implements `num::Float`, and supports all ordinary floating point operations.
/// 
/// This allows simple implementation of special floating point types, such as TensorFloat, IEEE754 Quadruple/binary128, Fp80, and BFloat16.
/// 
/// The accuracy of all of the floating point operations are not perfect, but work well enough to be usable. Various plots showing the accuracy of basic functions are shown in the [plots](https://github.com/sigurd4/custom_float/tree/master/plots) subfolder.
/// 
/// All floats can be converted into each other painlessly, though the conversion may produce rounding errors or unbounded outputs when converting to a float with lesser resolution.
///
/// # Examples
///
/// ```rust
/// #![feature(generic_const_exprs)]
///
/// use custom_float::Fp;
///
/// type FpSingle = Fp<u32, true, 8, 0, 23, 2>;
///
/// let two = FpSingle::from(2);
/// let four = FpSingle::from(4);
/// 
/// assert_eq!(two + two, four);
/// ```
#[derive(Clone, Copy)]
pub struct Fp<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize>(U)
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:;

impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:
{
    /// Size of floating-point number in bits
    pub const BIT_SIZE: usize = Self::SIGN_SIZE + Self::EXP_SIZE + Self::INT_SIZE + Self::FRAC_SIZE;
    /// Size of the sign bit
    pub const SIGN_SIZE: usize = Self::SIGN_BIT as usize;

    /// Wether or not the number has a sign bit. If not, it can only be positive.
    pub const SIGN_BIT: bool = SIGN_BIT;
    /// The size of the exponent part in bits
    pub const EXP_SIZE: usize = EXP_SIZE;
    /// The size of the integer part in bits
    pub const INT_SIZE: usize = INT_SIZE;
    /// The size of the fractional part in bits
    pub const FRAC_SIZE: usize = FRAC_SIZE;
    /// The base for the exponent
    pub const EXP_BASE: usize = EXP_BASE;

    /// Position of the sign bit
    pub const SIGN_POS: usize = Self::EXP_SIZE + Self::INT_SIZE + Self::FRAC_SIZE;
    /// Position of the first exponent bit
    pub const EXP_POS: usize = Self::INT_SIZE + Self::FRAC_SIZE;
    /// Position of the first integer bit
    pub const INT_POS: usize = Self::FRAC_SIZE;
    /// Position of the first fractional bit
    pub const FRAC_POS: usize = 0;

    /// Number of significant digits in base 2.
    pub const MANTISSA_DIGITS: usize = Self::INT_SIZE + Self::FRAC_SIZE;

    /// `true` if the number contains an implicit integer bit
    pub const IS_INT_IMPLICIT: bool = Self::INT_SIZE == 0;

    const MANTISSA_OP_SIZE: usize = FRAC_SIZE + INT_SIZE + Self::IS_INT_IMPLICIT as usize;
    const BASE_PADDING: usize = util::bitsize_of::<usize>() - EXP_BASE.leading_zeros() as usize - 1;

    #[cfg(any(test, feature = "use_std_float"))]
    fn as_lossless<const N: usize, const M: usize>(
        value: [Self; N],
        as_f16: impl FnOnce([f16; N]) -> [f16; M],
        as_f32: impl FnOnce([f32; N]) -> [f32; M],
        as_f64: impl FnOnce([f64; N]) -> [f64; M],
        as_f128: impl FnOnce([f128; N]) -> [f128; M]
    ) -> Option<[Self; M]>
    {
        macro_rules! as_float {
            ($($t:ty => $fn:expr),*) => {
                $(
                    if util::is_float_conversion_mutually_lossless::<Self, $t>()
                    {
                        return Some($fn(value.map(<$t as From<Self>>::from)).map(Self::from))
                    }
                )*
            };
        }

        as_float!(
            f16 => as_f16,
            f32 => as_f32,
            f64 => as_f64,
            f128 => as_f128
        );

        None
    }

    pub fn approx_eq(mut self, mut rhs: Self) -> bool
    {
        if self == rhs
            || (self.is_nan() && rhs.is_nan())
            || ((self.is_zero() || self.is_subnormal()) && (rhs.is_zero() || rhs.is_subnormal()))
            || (self.is_infinite() && rhs.is_infinite() && self.is_sign_negative() == rhs.is_sign_negative())
        {
            return true
        }

        if let Some(mut e0) = self.exp_bits().to_usize() && let Some(f0) = <u128 as NumCast>::from(self.mantissa_bits())
            && let Some(mut e1) = rhs.exp_bits().to_usize() && let Some(f1) = <u128 as NumCast>::from(rhs.mantissa_bits())
        {
            if !Self::IS_INT_IMPLICIT
            {
                self.normalize_up();
                rhs.normalize_up();
            }

            let tol = Self::FRAC_SIZE.saturating_sub(4);

            if self.is_sign_negative() == rhs.is_sign_negative()
            {
                if EXP_BASE != 2
                {
                    e0 = util::count_digits_in_base(e0, EXP_BASE);
                    e1 = util::count_digits_in_base(e1, EXP_BASE);
                }
                let a = (e0 + f0.ilog2() as usize)/tol;
                let b = (e1 + f1.ilog2() as usize)/tol;
    
                if a == b
                {
                    return true
                }
            }
            else
            {
                let mut e = e0 + e1;
                if EXP_BASE != 2
                {
                    e = util::count_digits_in_base(e, EXP_BASE);
                }
                let a = (e + (f0 + f1).ilog2() as usize)/tol;
    
                if a == 0
                {
                    return true
                }
            }
        }

        false
    }

    #[inline]
    fn max_exponent_bits() -> U
    {
        match EXP_SIZE < util::bitsize_of::<U>()
        {
            true => (U::one() << EXP_SIZE) - U::one(),
            false => U::max_value()
        }
    }
    #[inline]
    fn max_frac_bits() -> U
    {
        match FRAC_SIZE < util::bitsize_of::<U>()
        {
            true => (U::one() << FRAC_SIZE) - U::one(),
            false => U::max_value()
        }
    }
    #[inline]
    fn max_mantissa_bits() -> U
    {
        match Self::MANTISSA_DIGITS < util::bitsize_of::<U>()
        {
            true => (U::one() << Self::MANTISSA_DIGITS) - U::one(),
            false => U::max_value()
        }
    }

    fn shift_sign(s: U) -> U
    {
        match SIGN_BIT
        {
            true => s << Self::SIGN_POS,
            false => Zero::zero()
        }
    }
    fn shift_exp(e: U) -> U
    {
        match EXP_SIZE != 0
        {
            true => e << Self::EXP_POS,
            false => Zero::zero()
        }
    }
    fn shift_int(e: U) -> U
    {
        match INT_SIZE != 0
        {
            true => e << Self::INT_POS,
            false => Zero::zero()
        }
    }
    fn shift_frac(f: U) -> U
    {
        match FRAC_SIZE != 0
        {
            true => f << Self::FRAC_POS,
            false => Zero::zero()
        }
    }

    fn mantissa_bits(&self) -> U
    {
        let f = self.frac_bits();
        if self.is_subnormal()
        {
            f << 1usize
        }
        else
        {
            f + (self.int_bits() << Self::INT_POS)
        }
    }

    fn mantissa_bits_abnormal(&self) -> U
    {
        let f = self.frac_bits();
        if self.is_subnormal()
        {
            f
        }
        else
        {
            f + (self.int_bits() << Self::INT_POS)
        }
    }

    fn unexplicit_int(mantissa: &mut U)
    {
        if Self::IS_INT_IMPLICIT
        {
            let one = U::one() << Self::INT_POS;
            *mantissa = *mantissa - one;
            debug_assert!(*mantissa < one)
        }
    }

    fn from_sign_exp_mantissa(sign: bool, exp: U, mantissa: U) -> Self
    {
        if !SIGN_BIT && sign
        {
            if mantissa.is_zero()
            {
                return Self::zero()
            }
            return Self::qnan()
        }
        Self::from_exp_mantissa(exp, mantissa).with_sign(sign)
    }

    fn from_exp_mantissa(exp: U, mut mantissa: U) -> Self
    {
        if mantissa.is_zero()
        {
            return Self::zero()
        }
        let bits = if exp.is_zero() && Self::IS_INT_IMPLICIT // subnormal
        {
            Self::shift_frac(util::rounding_div_2(mantissa))
        }
        else
        {
            if exp >= Self::max_exponent_bits()
            {
                return Self::infinity()
            }

            Self::unexplicit_int(&mut mantissa);
            Self::shift_exp(exp) + Self::shift_frac(mantissa)
        };
        Self::from_bits(bits)
    }

    fn abs_add_mantissas(exp: &mut U, mut mantissa1: U, mut mantissa2: U, neg: bool) -> U
    {
        if neg
        {
            if mantissa1 >= mantissa2
            {
                mantissa1 - mantissa2
            }
            else
            {
                mantissa2 - mantissa1
            }
        }
        else 
        {
            match mantissa1.checked_add(&mantissa2)
            {
                Some(f) => f,
                None => {
                    if let Some(base) = U::from(EXP_BASE)
                    {
                        *exp = *exp + U::one();
                        mantissa1 = util::rounding_div(mantissa1, base);
                        mantissa2 = util::rounding_div(mantissa2, base);
                        mantissa1 + mantissa2
                    }
                    else if let Some(base) = U::from(EXP_BASE/2)
                    {
                        *exp = *exp + U::one();
                        mantissa1 = util::rounding_div(mantissa1, base);
                        mantissa2 = util::rounding_div(mantissa2, base);
                        util::rounding_div(mantissa1 + mantissa2, base)
                    }
                    else
                    {
                        U::max_value()
                    }
                }
            }
        }
    }

    fn add_signs(sign1: bool, sign2: bool, mantissa1: U, mantissa2: U) -> bool
    {
        match (sign1, sign2)
        {
            (false, false) => false,
            (false, true) => mantissa1 < mantissa2,
            (true, false) => mantissa1 > mantissa2,
            (true, true) => true,
        }
    }

    fn max_exponents(exp1: U, exp2: U, mantissa1: &mut U, mantissa2: &mut U) -> U
    {
        match exp1.cmp(&exp2)
        {
            Ordering::Less => {
                let shr = exp2 - exp1;
                *mantissa1 = if let Some(base) = U::from(EXP_BASE)
                {
                    util::rounding_div_pow(*mantissa1, base, shr)
                }
                else
                {
                    U::zero()
                };
                exp2
            },
            Ordering::Equal => exp1,
            Ordering::Greater => {
                let shr = exp1 - exp2;
                *mantissa2 = if let Some(base) = U::from(EXP_BASE)
                {
                    util::rounding_div_pow(*mantissa2, base, shr)
                }
                else
                {
                    U::zero()
                };
                exp1
            }
        }
    }

    fn add_nan(self, rhs: Self) -> Self
    {
        if rhs.is_snan() || !self.is_nan()
        {
            rhs
        }
        else
        {
            self
        }
    }

    fn add_with_sign(self, rhs: Self, neg: bool) -> Self
    {
        as_lossless!(
            [self, rhs],
            |[lhs, rhs]| [if neg {lhs - rhs} else {lhs + rhs}],
            {
                if self.is_nan() || rhs.is_nan()
                {
                    return self.add_nan(rhs)
                }
                
                if rhs.is_zero()
                {
                    return self
                }
                if self.is_zero()
                {
                    return if neg {-rhs} else {rhs}
                }

                let s0 = self.is_sign_negative();
                let s1 = rhs.is_sign_negative() ^ neg;

                match (self.is_infinite(), rhs.is_infinite())
                {
                    (true, true) => return if s0 == s1
                    {
                        self
                    }
                    else
                    {
                        Self::qnan()
                    },
                    (true, false) => return self,
                    (false, true) => return if neg {-rhs} else {rhs},
                    (false, false) => ()
                }
                
                let e0 = self.exp_bits();
                let e1 = rhs.exp_bits();
                let mut f0 = self.mantissa_bits();
                let mut f1 = rhs.mantissa_bits();

                let mut e = Self::max_exponents(e0, e1, &mut f0, &mut f1);
                let s = Self::add_signs(s0, s1, f0, f1);
                let mut f = Self::abs_add_mantissas(&mut e, f0, f1, s0 != s1);

                if f.is_zero()
                {
                    return if s {-Self::zero()} else {Self::zero()}
                }

                Self::normalize_mantissa(&mut e, &mut f, None);
                Self::from_sign_exp_mantissa(s, e, f)
            }
        )
    }

    /// Converts to this type from the input type.
    #[inline]
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn from<T>(from: T) -> Self
    where
        Self: From<T>
    {
        <Self as From<T>>::from(from)
    }

    fn normalize_up(&mut self)
    {
        let s = self.is_sign_negative();
        let mut e = self.exp_bits();
        let mut m = self.mantissa_bits();
        Self::normalize_mantissa_up(&mut e, &mut m, None);
        *self = Self::from_exp_mantissa(e, m).with_sign(s)
    }
    /*fn normalize_down(&mut self)
    {
        let s = self.is_sign_negative();
        let mut e = self.exp_bits();
        let mut m = self.mantissa_bits();
        Self::normalize_mantissa_down(&mut e, &mut m, None);
        *self = Self::from_exp_mantissa(e, m).with_sign(s)
    }
    fn normalize(&mut self)
    {
        let s = self.is_sign_negative();
        let mut e = self.exp_bits();
        let mut m = self.mantissa_bits();
        Self::normalize_mantissa(&mut e, &mut m, None);
        *self = Self::from_exp_mantissa(e, m).with_sign(s)
    }*/

    fn normalize_mantissa_down<M: AnyInt>(exp: &mut U, mantissa: &mut M, target_exp: Option<U>)
    {
        if mantissa.is_zero()
        {
            *exp = U::zero();
            return
        }
        let base = M::from(EXP_BASE);
        if let Some(base) = base
        {
            let min_exp = target_exp.unwrap_or_else(U::zero);
            while *exp > min_exp && (*mantissa >> (Self::MANTISSA_OP_SIZE - Self::BASE_PADDING)).is_zero()
            {
                *exp = *exp - U::one();
                *mantissa = *mantissa*base;
            }
        }
    }
    fn normalize_mantissa_up<M: AnyInt>(exp: &mut U, mantissa: &mut M, target_exp: Option<U>)
    {
        if mantissa.is_zero()
        {
            *exp = U::zero();
            return
        }
        if let Some(base) = M::from(EXP_BASE)
        {
            let max_exp = target_exp.unwrap_or_else(Self::max_exponent_bits);
            while *exp < max_exp && !(*mantissa >> Self::MANTISSA_OP_SIZE).is_zero()
            {
                *exp = *exp + U::one();
                *mantissa = util::rounding_div(*mantissa, base)
            }
        }
        else
        {
            *exp = U::zero();
            *mantissa = M::zero()
        }
    }

    fn normalize_mantissa<M: AnyInt>(exp: &mut U, mantissa: &mut M, target_exp: Option<U>)
    {
        Self::normalize_mantissa_down(exp, mantissa, target_exp);
        Self::normalize_mantissa_up(exp, mantissa, target_exp);
    }

    /// Converts from one custom floating-point number to another.
    /// Rounding errors may occurr.
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn from_fp<V: UInt, const S: bool, const E: usize, const I: usize, const F: usize, const B: usize>(fp: Fp<V, S, E, I, F, B>) -> Self
    where
        [(); util::bitsize_of::<V>() - S as usize - E - I - F]:,
        [(); B - 2]:
    {
        let s = fp.is_sign_negative();

        if S && !SIGN_BIT && s
        {
            return Self::qnan()
        }

        if EXP_SIZE == E && INT_SIZE == I && EXP_BASE == B && (I == 0 || FRAC_SIZE >= F)
        {
            if let Some(b) = if SIGN_BIT == S && FRAC_SIZE == F
            {
                <U as NumCast>::from(fp.to_bits())
            }
            else if util::bitsize_of::<U>() >= util::bitsize_of::<V>()
            {
                <U as NumCast>::from(fp.to_bits())
                    .map(|b| if FRAC_SIZE >= F
                    {
                        b << (FRAC_SIZE - F)
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
                    fp.to_bits() << (FRAC_SIZE - F)
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

        match fp.classify()
        {
            FpCategory::Nan => {
                if fp.is_snan()
                {
                    Self::snan()
                }
                else
                {
                    Self::qnan()
                }
            },
            FpCategory::Infinite => Self::infinity(),
            FpCategory::Zero => Self::zero(),
            FpCategory::Subnormal | FpCategory::Normal => {
                let mut e1 = fp.exp_bits();
                let mut f = fp.mantissa_bits();

                let df = FRAC_SIZE as isize - F as isize;

                let base1 = V::from(B);

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
                            f = if let Some(base1) = base1
                            {
                                util::rounding_div(f, base1)
                            }
                            else
                            {
                                V::zero()
                            }
                        }
                    }
                };
                
                let bias1 = Fp::<V, S, E, I, F, B>::exp_bias();
                let bias2 = Self::exp_bias();

                let mut e = {
                    match if EXP_BASE == B
                    {
                        if util::bitsize_of::<U>() < util::bitsize_of::<V>() && let Some(bias) = V::from(bias2)
                        {
                            match bias1.cmp(&bias)
                            {
                                Ordering::Greater => e1.checked_sub(&(bias1 - bias)),
                                Ordering::Equal => Some(e1),
                                Ordering::Less => e1.checked_add(&(bias - bias1))
                            }.and_then(U::from)
                        }
                        else if let Some(bias) = U::from(bias1) && let Some(e) = U::from(e1)
                        {
                            match bias.cmp(&bias2)
                            {
                                Ordering::Greater => e.checked_sub(&(bias - bias2)),
                                Ordering::Equal => Some(e),
                                Ordering::Less => e.checked_add(&(bias2 - bias))
                            }
                        }
                        else
                        {
                            None
                        }
                    }
                    else
                    {
                        None
                    }
                    {
                        Some(e) => e,
                        None => {
                            let base1 = U::from(B);
                            let base2 = U::from(EXP_BASE);
        
                            let mut e = bias2;
                            if let Some(base1) = base1
                            {
                                while e1 > bias1
                                {
                                    e1 = e1 - V::one();
                                    while f.leading_zeros() as usize <= Fp::<V, S, E, I, F, B>::BASE_PADDING
                                    {
                                        e = e + U::one();
                                        f = if let Some(base2) = base2
                                        {
                                            util::rounding_div(f, base2)
                                        }
                                        else
                                        {
                                            U::zero()
                                        }
                                    }
                                    f = f*base1
                                }
                            }
                            if let Some(base2) = base2
                            {
                                while e1 < bias1
                                {
                                    e1 = e1 + V::one();
                                    while e > U::zero() && f.leading_zeros() as usize > Self::BASE_PADDING
                                    {
                                        e = e - U::one();
                                        f = f*base2;
                                    }
                                    f = if let Some(base1) = base1
                                    {
                                        util::rounding_div(f, base1)
                                    }
                                    else
                                    {
                                        U::zero()
                                    };
                                }
                            }
                            e
                        }
                    }
                };

                Self::normalize_mantissa(&mut e, &mut f, None);
                Self::from_exp_mantissa(e, f)
            },
        }.with_sign(s)
    }

    /// Converts an unsigned integer into a custom floating-point type.
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn from_uint<I: UInt>(from: I) -> Self
    {
        let mut e = Self::exp_bias();
        
        let f = if util::bitsize_of::<I>() > util::bitsize_of::<U>()
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

            Self::normalize_mantissa(&mut e, &mut f, None);
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
            
            Self::normalize_mantissa(&mut e, &mut f, None);
            f
        };

        Self::from_exp_mantissa(e, f)
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
        if s && !SIGN_BIT
        {
            return Self::qnan()
        }
        let mut e = Self::exp_bias();
        let f = if util::bitsize_of::<I>() - 1 > util::bitsize_of::<U>()
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

            Self::normalize_mantissa(&mut e, &mut f, None);
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

            Self::normalize_mantissa(&mut e, &mut f, None);
            f
        };

        Self::from_sign_exp_mantissa(s, e, f)
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
        let f = self.mantissa_bits();
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

            I::from(f)?
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
        let f = self.mantissa_bits();
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
        let f = self.mantissa_bits();
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

            I::from(f)?
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
        let f = self.mantissa_bits();
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
            else if f >= U::from(I::max_value()).unwrap()
            {
                return I::max_value()
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
    #[inline]
    pub fn sign_bit(self) -> U
    {
        if !SIGN_BIT
        {
            return U::zero()
        }
        (self.to_bits() & U::max_value() >> (util::bitsize_of::<U>() - Self::SIGN_POS - Self::SIGN_SIZE)) >> Self::SIGN_POS
    }
    /// Returns the exponent bits of the custom floating-point number.
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub fn exp_bits(self) -> U
    {
        if EXP_SIZE == 0
        {
            return U::zero()
        }
        (self.to_bits() & U::max_value() >> (util::bitsize_of::<U>() - Self::EXP_POS - EXP_SIZE)) >> Self::EXP_POS
    }
    /// Returns the integer bits of the custom floating-point number.
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub fn int_bits(self) -> U
    {
        if Self::IS_INT_IMPLICIT
        {
            return if self.is_normal() {U::one()} else {U::zero()}
        }
        (self.to_bits() & U::max_value() >> (util::bitsize_of::<U>() - Self::INT_POS - INT_SIZE)) >> Self::INT_POS
    }
    /// Returns the fractional bits of the custom floating-point number.
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub fn frac_bits(self) -> U
    {
        if FRAC_SIZE == 0
        {
            return U::zero()
        }
        (self.to_bits() & U::max_value() >> (util::bitsize_of::<U>() - Self::FRAC_POS - FRAC_SIZE)) >> Self::FRAC_POS
    }

    /// Returns the exponent bias
    #[must_use]
    #[inline]
    pub fn exp_bias() -> U
    {
        if EXP_SIZE == 0
        {
            return U::zero()
        }
        U::max_value() >> (util::bitsize_of::<U>() + 1 - EXP_SIZE)
    }
    
    /// Returns the `NaN` value.
    ///
    /// # Examples
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let nan = FpSingle::nan();
    ///
    /// assert!(nan.is_nan());
    /// ```
    #[must_use]
    #[inline]
    pub fn nan() -> Self
    {
        Self::snan()
    }
    
    /// Returns the `qNaN` value.
    ///
    /// # Examples
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let qnan = FpSingle::qnan();
    ///
    /// assert!(qnan.is_nan());
    /// assert!(!qnan.is_snan());
    /// ```
    #[must_use]
    #[inline]
    pub fn qnan() -> Self
    {
        Self::from_bits(U::max_value() >> (util::bitsize_of::<U>() - Self::SIGN_POS))
    }
    
    /// Returns the `sNaN` value.
    ///
    /// # Examples
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let snan = FpSingle::snan();
    ///
    /// assert!(snan.is_nan());
    /// assert!(snan.is_snan());
    /// ```
    #[must_use]
    #[inline]
    pub fn snan() -> Self
    {
        if Self::INT_POS + INT_SIZE < 1
        {
            return Self::qnan()
        }
        Self::from_bits((U::max_value() >> (util::bitsize_of::<U>() - Self::SIGN_POS)) - (U::one() << (Self::INT_POS + INT_SIZE - 1)))
    }

    /// Returns `true` if the number is a signaling NaN.
    ///
    /// # Examples
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let snan = FpSingle::snan();
    /// let qnan = FpSingle::qnan();
    ///
    /// assert!(snan.is_snan());
    /// assert!(!qnan.is_snan());
    /// ```
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub fn is_snan(self) -> bool
    {
        if Self::INT_POS + INT_SIZE < 1 || !self.is_nan()
        {
            return false
        }
        (self.to_bits() >> (Self::INT_POS + INT_SIZE - 1)) & U::one() == U::zero()
    }

    /// Returns the infinite value.
    ///
    /// # Examples
    ///
    /// ```
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let infinity = FpSingle::infinity();
    ///
    /// assert!(infinity.is_infinite());
    /// assert!(!infinity.is_finite());
    /// assert!(infinity > FpSingle::max_value());
    /// ```
    #[must_use]
    #[inline]
    pub fn infinity() -> Self
    {
        Self::from_bits(Self::shift_exp(U::max_value() >> (util::bitsize_of::<U>() - EXP_SIZE)))
    }

    /// Returns the negative infinite value.
    ///
    /// # Examples
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
    #[inline]
    pub fn neg_infinity() -> Self
    {
        -Self::infinity()
    }

    /// Returns `-0.0`.
    ///
    /// # Examples
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
    #[inline]
    pub fn neg_zero() -> Self
    {
        -Self::zero()
    }

    /// Returns the smallest finite value.
    ///
    /// # Examples
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
    #[inline]
    pub fn min_value() -> Self
    {
        if !SIGN_BIT
        {
            return Self::zero()
        }
        -Self::max_value()
    }

    /// Returns the smallest positive, normal value.
    ///
    /// # Examples
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
    #[inline]
    pub fn min_positive_value() -> Self
    {
        let mut bits = Self::shift_exp(U::one());
        if !Self::IS_INT_IMPLICIT
        {
            bits = bits + Self::shift_int(U::one());
        }
        Self::from_bits(bits)
    }

    /// [Machine epsilon] value.
    ///
    /// This is the difference between `1.0` and the next larger representable number.
    ///
    /// [Machine epsilon]: https://en.wikipedia.org/wiki/Machine_epsilon
    ///
    /// # Examples
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
    #[inline]
    pub fn epsilon() -> Self
    {
        let bias = Self::exp_bias();
        
        if !Self::IS_INT_IMPLICIT
        {
            return Self::from_bits(Self::shift_exp(bias) + Self::shift_frac(U::one()))
        }

        let exp_frac = U::from(util::count_digits_in_base(FRAC_SIZE + INT_SIZE, EXP_BASE)).unwrap();
        if bias <= exp_frac
        {
            return Self::from_bits(util::powu(U::from(EXP_BASE).unwrap(), <usize as NumCast>::from(bias).unwrap() - 1))
        }
        Self::from_bits(Self::shift_exp(bias - exp_frac))
    }

    /// Returns the largest finite value.
    ///
    /// # Examples
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
    #[inline]
    pub fn max_value() -> Self
    {
        Self::from_bits(Self::shift_exp(Self::max_exponent_bits()) - U::one())
    }

    /// Returns `true` if this value is NaN.
    ///
    /// # Examples
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
    #[inline]
    pub fn is_nan(self) -> bool
    {
        !self.is_finite() && !(self.frac_bits().is_zero() && (INT_SIZE == 0 || self.int_bits().is_zero()))
    }

    /// Returns `true` if this value is positive infinity or negative infinity, and
    /// `false` otherwise.
    ///
    /// # Examples
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
    #[inline]
    pub fn is_infinite(self) -> bool
    {
        !self.is_finite() && self.frac_bits().is_zero() && (INT_SIZE == 0 || self.int_bits().is_zero())
    }

    /// Returns `true` if this number is neither infinite nor NaN.
    ///
    /// # Examples
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
    #[inline]
    pub fn is_finite(self) -> bool
    {
        self.exp_bits() != Self::max_exponent_bits()
    }

    /// Returns `true` if the number is neither zero, infinite,
    /// [subnormal], or NaN.
    ///
    /// # Examples
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
    #[inline]
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
    /// # Examples
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
    #[inline]
    pub fn is_subnormal(self) -> bool
    {
        Self::IS_INT_IMPLICIT && self.exp_bits() == U::zero() && !self.is_zero()
    }

    /// Returns the floating point category of the number. If only one property
    /// is going to be tested, it is generally faster to use the specific
    /// predicate instead.
    ///
    /// # Examples
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
    #[inline]
    pub fn classify(self) -> FpCategory
    {
        let e = self.exp_bits();
        if self.is_zero()
        {
            return FpCategory::Zero
        }
        if Self::IS_INT_IMPLICIT && e.is_zero()
        {
            return FpCategory::Subnormal
        }
        if e == (U::max_value() >> (util::bitsize_of::<U>() - EXP_SIZE))
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
    /// # Examples
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
    #[inline]
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
    /// # Examples
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
    #[inline]
    pub fn ceil(mut self) -> Self
    {
        if !self.is_finite()
        {
            return self
        }
        let mut m = self % Self::one();
        if self.is_sign_positive() && !m.is_zero()
        {
            if !SIGN_BIT
            {
                self += Self::one()
            }
            else
            {
                m -= Self::one()
            }
        }
        self - m
    }

    /// Returns the nearest integer to `self`. If a value is half-way between two
    /// integers, round away from `0.0`.
    ///
    /// # Examples
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
    #[inline]
    pub fn round(mut self) -> Self
    {
        if !self.is_finite()
        {
            return self
        }
        let mut m = self % Self::one();
        let half = Self::from(0.5);
        if self.is_sign_positive() && m >= half
        {
            if !SIGN_BIT
            {
                self += Self::one()
            }
            else
            {
                m -= Self::one()
            }
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
    #[inline]
    pub fn round_ties_even(mut self) -> Self
    {
        if !self.is_finite()
        {
            return self
        }
        let one = Self::one();
        let t = self % Self::from(2.0);
        let mut m = self % one;
        let half = Self::from(0.5);
        if self.is_sign_positive() && m >= half && (m != half || t > one)
        {
            if !SIGN_BIT
            {
                self += Self::one()
            }
            else
            {
                m -= Self::one()
            }
        }
        if self.is_sign_negative() && m <= -half && (m != -half || t < -one)
        {
            m += Self::one()
        }
        self - m
    }

    /// Returns the integer part of `self`.
    /// This means that non-integer numbers are always truncated towards zero.
    ///
    /// # Examples
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
    #[inline]
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
    /// # Examples
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
    #[inline]
    pub fn fract(self) -> Self
    {
        self - self.trunc()
    }

    /// Computes the absolute value of `self`. Returns `NaN` if the number is `NaN`.
    ///
    /// # Examples
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
    #[inline]
    pub fn abs(self) -> Self
    {
        if !SIGN_BIT
        {
            return self
        }
        Self::from_bits(self.to_bits() & (U::max_value() >> (util::bitsize_of::<U>() - Self::SIGN_POS)))
    }

    /// Returns a number that represents the sign of `self`.
    ///
    /// - `1.0` if the number is positive, `+0.0` or `inf`
    /// - `-1.0` if the number is negative, `-0.0` or `-inf`
    /// - `NaN` if the number is `NaN`
    ///
    /// # Examples
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
    #[inline]
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
    /// # Examples
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
    #[inline]
    pub fn is_sign_positive(self) -> bool
    {
        if !SIGN_BIT
        {
            return true
        }
        self.sign_bit().is_zero()
    }

    /// Returns `true` if `self` has a negative sign, including `-0.0`, NaNs with
    /// negative sign bit and negative infinity. Note that IEEE 754 doesn't assign any
    /// meaning to the sign bit in case of a NaN, and as Rust doesn't guarantee that
    /// the bit pattern of NaNs are conserved over arithmetic operations, the result of
    /// `is_sign_negative` on a NaN might produce an unexpected result in some cases.
    /// See [explanation of NaN as a special value](f32) for more info.
    ///
    /// # Examples
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
    #[inline]
    pub fn is_sign_negative(self) -> bool
    {
        if !SIGN_BIT
        {
            return false
        }
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
    /// # Examples
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
        self.next(false)
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
    /// # Examples
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
        self.next(true)
    }

    fn next(self, sign: bool) -> Self
    {
        if self.is_nan()
        {
            return self
        }

        let mut s = self.is_sign_negative();

        if self.is_infinite()
        {
            if s ^ sign
            {
                return Self::max_value().with_sign(!sign)
            }
            return self
        }

        let mut e = self.exp_bits();

        if Self::MANTISSA_DIGITS == 0
        {
            if s ^ sign
            {
                if e.is_zero()
                {
                    return -self
                }
                else
                {
                    e = e - U::one()
                }
            }
            else
            {
                e = e + U::one()
            };
            let s_bits = if s {Self::shift_sign(U::one())} else {U::zero()};
            Self::from_bits(s_bits + Self::shift_exp(e))
        }
        else if Self::IS_INT_IMPLICIT
        {
            let mut f = self.frac_bits();

            if s ^ sign
            {
                if f.is_zero()
                {
                    if e.is_zero()
                    {
                        s = !s;
                    }
                    else
                    {
                        e = e - U::one();
                        f = Self::max_frac_bits()
                    }
                }
                else
                {
                    f = f - U::one()
                }
            }
            else if f == Self::max_frac_bits()
            {
                e = e + U::one();
                f = U::zero()
            }
            else
            {
                f = f + U::one()
            }

            let s_bits = if s {Self::shift_sign(U::one())} else {U::zero()};
            Self::from_bits(s_bits + Self::shift_exp(e) + Self::shift_frac(f))
        }
        else
        {
            let mut f = self.mantissa_bits_abnormal();

            let base = U::from(EXP_BASE);

            if s ^ sign
            {
                if f.is_zero()
                {
                    s = !s;
                    e = U::zero()
                }
                else if let Some(base) = base && let Some(mut ff) = f.checked_mul(&base)
                {
                    let mut o = U::one();
                    while !e.is_zero() && let oo = ff.saturating_sub(Self::max_mantissa_bits()) && oo < base
                    {
                        o = oo.max(U::one());
                        e = e - U::one();
                        f = ff;
                        if let Some(fff) = f.checked_mul(&base)
                        {
                            ff = fff
                        }
                        else
                        {
                            break
                        }
                    }
                    f = f - o
                }
                else if !e.is_zero() && f == U::one()
                {
                    e = e - U::one();
                    f = Self::max_mantissa_bits()
                }
                else
                {
                    f = f - U::one()
                }
            }
            else
            {
                if let Some(base) = base && let Some(mut ff) = f.checked_mul(&base)
                {
                    while !e.is_zero() && ff <= Self::max_mantissa_bits()
                    {
                        e = e - U::one();
                        f = ff;
                        if let Some(fff) = f.checked_mul(&base)
                        {
                            ff = fff
                        }
                        else
                        {
                            break
                        }
                    }
                }
                if f == Self::max_mantissa_bits()
                {
                    e = e + U::one();
                    if e == Self::max_exponent_bits()
                    {
                        let s_bits = if s {Self::shift_sign(U::one())} else {U::zero()};
                        return Self::from_bits(s_bits + Self::shift_exp(e))
                    }
                    if let Some(base) = base
                    {
                        f = f/base + U::one();
                    }
                    else
                    {
                        f = U::one()
                    }
                }
                else
                {
                    f = f + U::one()
                }
            }

            let s_bits = if s {Self::shift_sign(U::one())} else {U::zero()};
            Self::from_bits(s_bits + Self::shift_exp(e) + Self::shift_frac(f))
        }
    }
    
    /// Returns the maximum of the two numbers, propagating NaN.
    ///
    /// This returns NaN when *either* argument is NaN, as opposed to
    /// [`Fp::max`] which only returns NaN when *both* arguments are NaN.
    ///
    /// # Examples
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
    /// # Examples
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
    #[inline]
    pub fn midpoint(self, other: Self) -> Self
    {
        let half = Self::from(0.5);

        let lo = Self::min_positive_value()*Self::from_uint(2u8);
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
    /// # Examples
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
    #[inline]
    pub fn mul_add(self, a: Self, b: Self) -> Self
    {
        (self*a) + b
    }

    fn integral_div(mut mantissa1: U, mut mantissa2: U, exp_offset: &mut U) -> Result<U, Self>
    {
        let base = U::from(EXP_BASE);
        loop
        {
            if mantissa2.is_zero()
            {
                return Err(Self::infinity())
            }
            let f = util::rounding_div(mantissa1, mantissa2);
            if let Some(base) = base && !mantissa1.is_zero() && f.leading_zeros() as usize >= Self::BASE_PADDING
            {
                if mantissa1.leading_zeros() as usize > Self::BASE_PADDING
                {
                    mantissa1 = mantissa1*base
                }
                else
                {
                    if mantissa2 % base != U::zero()
                    {
                        break Ok(f)
                    }
                    mantissa2 = mantissa2/base;
                }
                *exp_offset = *exp_offset + U::one();
            }
            else
            {
                break Ok(f)
            }
        }
    }

    fn mantissa_div(mantissa1: U, mantissa2: U, exp: &mut U, exp_offset: &mut U) -> Result<U, Self>
    {
        let mut mantissa = Self::integral_div(mantissa1, mantissa2, exp_offset)?;
        Self::integral_to_mantissa(&mut mantissa, exp, exp_offset)?;
        Ok(mantissa)
    }

    fn exponent_sub(exp1: U, exp2: U, mantissa1: &mut U, mantissa2: &mut U) -> Result<U, Self>
    {
        let bias = Self::exp_bias();
        if let Some(e) = exp1.checked_sub(&exp2)
        {
            return e.checked_add(&bias).ok_or_else(Self::infinity)
        }
        if let Some(e) = bias.checked_sub(&exp2)
        {
            return exp1.checked_add(&e).ok_or_else(Self::infinity)
        }
        let e = exp2 - bias;
        if let Some(e) = exp1.checked_sub(&e)
        {
            return Ok(e)
        }
        if let Some(base) = U::from(EXP_BASE)
        {
            let mut o = e - exp1;
            while o > U::zero()
            {
                o = o - U::one();
                if (*mantissa1 % base).is_zero()
                {
                    *mantissa1 = *mantissa1/base;
                }
                else if mantissa2.leading_zeros() as usize > Self::BASE_PADDING
                {
                    *mantissa2 = *mantissa2*base;
                }
                else
                {
                    *mantissa1 = util::rounding_div(*mantissa1, base);
                }
            }
            return Ok(U::zero())
        }
        Err(Self::zero())
    }

    fn integral_to_mantissa(mantissa: &mut U, exp: &mut U, exp_offset: &mut U) -> Result<(), Self>
    {
        if FRAC_SIZE == 0
        {
            return Ok(())
        }
        if EXP_BASE.is_power_of_two() && let Some(change) = U::from(FRAC_SIZE/EXP_BASE.ilog2() as usize)
        {
            if let Some(diff) = exp_offset.checked_sub(&change)
            {
                *exp_offset = diff
            }
            else if let Some(diff) = exp.checked_add(&(change - *exp_offset))
            {
                *exp_offset = U::zero();
                *exp = diff
            }
            else
            {
                return Err(Self::infinity())
            }
        }
        else if let Some(base) = U::from(EXP_BASE)
        {
            for _ in 0..FRAC_SIZE
            {
                *mantissa = loop
                {
                    if mantissa.is_zero()
                    {
                        return Err(Self::zero())
                    }
                    if mantissa.leading_zeros() > 0
                    {
                        break *mantissa << 1usize
                    }
                    if *exp_offset > U::zero()
                    {
                        *exp_offset = *exp_offset - U::one();
                    }
                    else if *exp < U::max_value()
                    {
                        *exp = *exp + U::one();
                    }
                    else
                    {
                        return Err(Self::infinity())
                    }
                    *mantissa = util::rounding_div(*mantissa, base);
                }
            }
        }
        else
        {
            return Err(Self::zero())
        }
        Ok(())
    }

    /// Take the reciprocal (inverse) of a number, `1/x`.
    ///
    /// # Examples
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
    #[inline]
    pub fn recip(self) -> Self
    {
        as_lossless!(
            [self],
            |[x]| [x.recip()],
            {
                let s = self.is_sign_negative();
                match self.classify()
                {
                    FpCategory::Nan => self,
                    FpCategory::Zero => Self::infinity().with_sign(s),
                    FpCategory::Infinite => Self::zero().with_sign(s),
                    FpCategory::Normal | FpCategory::Subnormal => {
                        if self.abs().is_one()
                        {
                            return self
                        }
                
                        let mut e: U = Self::exp_bias() + Self::exp_bias() - self.exp_bits();
                        let f0: U = U::one() << Self::INT_POS;
                        let f1: U = self.mantissa_bits();
                        
                        let mut o = U::zero();
                        let mut f = match Self::mantissa_div(f0, f1, &mut e, &mut o)
                        {
                            Ok(f) => f,
                            Err(done) => return done.with_sign(s)
                        };
                        Self::normalize_mantissa_up(&mut e, &mut f, Some(o));
                        let mut e = match e.checked_sub(&o)
                        {
                            Some(e) => e,
                            None => return Self::zero().with_sign(s)
                        };
                
                        Self::normalize_mantissa(&mut e, &mut f, None);
                        Self::from_sign_exp_mantissa(s, e, f)
                    }
                }
            }
        )
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
    #[inline]
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
    #[inline]
    pub fn powu<I: UInt>(self, n: I) -> Self
    {
        util::powu(self, n)
    }

    /// Raises a number to a floating point power.
    ///
    /// This implementation is based on the [Apple Libm-315 implementation of powf](https://opensource.apple.com/source/Libm/Libm-315/Source/ARM/powf.c)
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
            let e = Self::max_exponent_bits() - U::one();
            let i = if INT_SIZE > 0
            {
                Self::shift_int(U::one())
            }
            else
            {
                U::zero()
            };
            Self::from_bits(Self::shift_exp(e) + i)
        };
        let edge_n = {
            let e = Self::exp_bias() + exp_frac;
            let i = if INT_SIZE > 0
            {
                Self::shift_int(U::one())
            }
            else
            {
                U::zero()
            };
            Self::from_bits(Self::shift_exp(e) + i)
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
        let n_d = (nabs.trunc() % Self::from(2u8)  + Self::one())*Self::from(0.5);
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
                return Self::neg_zero()
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
        if Self::from(0.5) == nabs
        {
            if self < Self::zero() || n > Self::zero()
            {
                return self.sqrt()
            }
            return (Self::one()/self).sqrt()
        }

        if !SIGN_BIT
        {
            let xabs_log = xabs.recip().logb();
    
            let n_xabs_log = n*xabs_log;
            
            n_xabs_log.expb().recip()
        }
        else
        {
            let xabs_log = xabs.logb();

            let n_xabs_log = n*xabs_log;
            
            n_xabs_log.expb()
        }
    }

    fn sqrt_y_naive(self) -> Self
    {
        if !SIGN_BIT && self < Self::one()
        {
            let xabs_log = self.recip().logb();

            let n_xabs_log = xabs_log/Self::from(2u8);
            
            n_xabs_log.expb().recip()
        }
        else
        {
            let xabs_log = self.logb();

            let n_xabs_log = xabs_log/Self::from(2u8);
            
            n_xabs_log.expb()
        }
    }

    /// Returns the square root of a number.
    ///
    /// Returns NaN if `self` is a negative number other than `-0.0`.
    ///
    /// This implementation is based on the fast sqrt described in: https://en.wikipedia.org/wiki/Methods_of_computing_square_roots#Approximations_that_depend_on_the_floating_point_representation
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
    #[inline]
    pub fn sqrt(self) -> Self
    {
        trait SqrtSpec: Sized
        {
            fn _y(self) -> Self; 
        }
        impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> SqrtSpec for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
        where
            [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
            [(); EXP_BASE - 2]:
        {
            default fn _y(self) -> Self
            {
                self.sqrt_y_naive()
            }
        }
        #[allow(clippy::identity_op)]
        impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> SqrtSpec for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
        where
            [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
            [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - 0 - FRAC_SIZE]:,
            [(); EXP_BASE - 2]:
        {
            fn _y(self) -> Self
            {
                if EXP_BASE != 2
                {
                    self.sqrt_y_naive()
                }
                else if !Self::IS_INT_IMPLICIT
                {
                    let x = Fp::<U, SIGN_BIT, EXP_SIZE, 0, FRAC_SIZE, EXP_BASE>::from_fp(self);
                    Self::from_fp(Fp::<U, SIGN_BIT, EXP_SIZE, 0, FRAC_SIZE, EXP_BASE>::from_bits(
                        ((Fp::<U, SIGN_BIT, EXP_SIZE, 0, FRAC_SIZE, EXP_BASE>::exp_bias() + U::one()) * (U::one() << (FRAC_SIZE - 1)))
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
                }
            }
        }

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

        let y = self._y();

        const NEWTON: usize = NEWTON_RT;
        let half = Self::from(0.5);
        let mut y = y;
        for _ in 0..NEWTON
        {
            let y_ = half*(y + self/y);
            if !y_.is_finite()
            {
                break
            }
            y = y_;
        }
        y
    }

    /// Returns `EXP_BASE^(self)`.
    ///
    /// This implementation is roughly based on the exp2 implementation described here: https://stackoverflow.com/questions/65554112/fast-double-exp2-function-in-c.
    ///
    /// # Examples
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
    /// let abs_difference_f = (f.expb() - FpDouble::from(4.0)).abs();
    ///
    /// // 10^2 - 100 == 0
    /// let abs_difference_d = (d.expb() - DecDouble::from(100.0)).abs();
    ///
    /// assert!(abs_difference_f < FpDouble::from(1e-10));
    /// assert!(abs_difference_d < DecDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn expb(self) -> Self
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
            return Self::from_bits((Self::max_exponent_bits()) << Self::EXP_POS)*self;
        }
        if self <= -Self::from_uint((U::one() << (EXP_SIZE - 1)) + exp_frac - U::one())
        {
            return -Self::from_bits(U::one() << (Self::EXP_POS - EXP_SIZE/2))*self;
        }

        let neg = self.is_sign_negative();
        let x = self.abs().max(-Self::from_uint(Self::exp_bias() - U::one()));
        
        let e = x.floor();
        let f = x - e;

        let ln_exp_base = if EXP_BASE == 2
        {
            Self::LN_2()
        }
        else if EXP_BASE == 10
        {
            Self::LN_10()
        }
        else
        {
            Self::from((EXP_BASE as f64).ln())
        };

        let z = Self::from(Into::<f64>::into(f*ln_exp_base).exp());

        let e = e + Self::from_uint(Self::exp_bias());
        if e > Self::from_uint(Self::max_exponent_bits())
        {
            return if neg {Self::zero()} else {Self::infinity()}*z
        }
        if e < Self::zero()
        {
            return if !neg {Self::zero()} else {Self::infinity()}*z
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
        let y = z*Self::from_bits(Self::shift_exp(e) + Self::shift_int(i));

        if neg
        {
            return y.recip()
        }
        y
    }

    #[inline]
    fn exp_nonewton(self) -> Self
    {
        if EXP_BASE == 2
        {
            return (self/Self::LN_2()).expb()
        }
        if EXP_BASE == 10
        {
            return (self/Self::LN_10()).expb()
        }
        
        (self/Self::from((EXP_BASE as f64).ln())).expb()
    }

    /// Returns `e^(self)`, (the exponential function).
    ///
    /// # Examples
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
    /// assert!(abs_difference < FpDouble::from(1e-9));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn exp(self) -> Self
    {
        let mut y = self.exp_nonewton();

        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_EXP;

            let one = Self::one();
            for _ in 0..NEWTON
            {
                let (x, x_s) = if SIGN_BIT || y >= one
                {
                    (y.ln_nonewton(), false)
                }
                else
                {
                    (y.recip().ln_nonewton(), true)
                };

                let (dxx, s) = util::add_extra_sign((x, x_s), (self, true), SIGN_BIT);
                let dy = y*dxx;

                if !dy.is_finite() || (!SIGN_BIT && !s && dy > y)
                {
                    break
                }
                if s
                {
                    y += dy
                }
                else
                {
                    y -= dy
                }
            }
        }

        y
    }

    /// Returns `10^(self)`.
    ///
    /// # Examples
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
    #[inline]
    pub fn exp10(self) -> Self
    {
        if EXP_BASE == 10
        {
            return self.expb()
        }
        (self*Self::LN_10()).exp()
    }

    /// Returns `2^(self)`.
    ///
    /// # Examples
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
    #[inline]
    pub fn exp2(self) -> Self
    {
        if EXP_BASE == 2
        {
            return self.expb()
        }
        (self*Self::LN_2()).exp()
    }

    #[inline]
    fn ln_nonewton(self) -> Self
    {
        if EXP_BASE == 2
        {
            return self.logb()*Self::LN_2()
        }
        if EXP_BASE == 10
        {
            return self.logb()*Self::LN_10()
        }
        self.logb()*Self::from((EXP_BASE as f64).ln())
    }

    /// Returns the natural logarithm of the number.
    ///
    /// # Examples
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
    /// assert!(abs_difference < FpDouble::from(1e-9));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn ln(self) -> Self
    {
        let mut y = self.ln_nonewton();

        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_LN;

            let one = Self::one();
            for _ in 0..NEWTON
            {
                let x = self/y.exp();
                let (dy, s) = if SIGN_BIT || one >= x
                {
                    (one - x, false)
                }
                else
                {
                    (x - one, true)
                };
                if !dy.is_finite() || (!SIGN_BIT && !s && dy > y)
                {
                    break
                }
                if s
                {
                    y += dy
                }
                else
                {
                    y -= dy
                }
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
    /// # Examples
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
    #[inline]
    pub fn log(self, base: Self) -> Self
    {
        self.ln()/base.ln()
    }

    /// Returns the logarithm base `EXP_BASE` of the number.
    ///
    /// # Examples
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
    /// let abs_difference_2 = (two.logb() - FpDouble::one()).abs();
    /// 
    /// // log10(10) - 1 == 0
    /// let abs_difference_10 = (ten.logb() - DecDouble::one()).abs();
    ///
    /// assert!(abs_difference_2 < FpDouble::from(1e-10));
    /// assert!(abs_difference_10 < DecDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn logb(self) -> Self
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
        
        let base = U::from(EXP_BASE).unwrap();
        if e.is_zero() && Self::IS_INT_IMPLICIT
        {
            while f < U::one() << (Self::MANTISSA_OP_SIZE - Self::BASE_PADDING)
            {
                y -= Self::one();
                f = f*base;
            }
        }

        Self::unexplicit_int(&mut f);

        let u = Self::from_bits(Self::shift_exp(bias) + Self::shift_frac(f));
        
        let u: f64 = u.into();
        y += Self::from(u.log(EXP_BASE as f64)); 
        y
    }
    
    #[inline]
    fn log2_nonewton(self) -> Self
    {
        if EXP_BASE == 2
        {
            return self.logb()
        }
        self.ln_nonewton()/Self::LN_2()
    }

    /// Returns the base 2 logarithm of the number.
    ///
    /// # Examples
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
    #[inline]
    pub fn log2(self) -> Self
    {
        let mut y = self.log2_nonewton();

        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_LN;

            let one = Self::one();
            let ln2 = Self::LN_2();
            for _ in 0..NEWTON
            {
                let xdyexp2 = self/y.exp2();
                let (dy, s) = if SIGN_BIT || one >= xdyexp2
                {
                    ((one - xdyexp2)/ln2, false)
                }
                else
                {
                    ((xdyexp2 - one)/ln2, true)
                };
                if !dy.is_finite() || (!SIGN_BIT && !s && dy > y)
                {
                    break
                }
                if s
                {
                    y += dy
                }
                else
                {
                    y -= dy
                }
            }
        }

        y
    }
    
    #[inline]
    fn log10_nonewton(self) -> Self
    {
        if EXP_BASE == 10
        {
            return self.logb()
        }
        self.ln_nonewton()/Self::LN_10()
    }

    /// Returns the base 10 logarithm of the number.
    ///
    /// # Examples
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
    /// assert!(abs_difference < FpDouble::from(1e-4));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn log10(self) -> Self
    {
        let mut y = self.log10_nonewton();

        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_LN;

            let one = Self::one();
            let ln10 = Self::LN_10();
            for _ in 0..NEWTON
            {
                let xdyexp10 = self/y.exp10();
                let (dy, s) = if SIGN_BIT || one >= xdyexp10
                {
                    ((one - xdyexp10)/ln10, false)
                }
                else
                {
                    ((xdyexp10 - one)/ln10, true)
                };
                if !dy.is_finite() || (!SIGN_BIT && !s && dy > y)
                {
                    break
                }
                if s
                {
                    y += dy
                }
                else
                {
                    y -= dy
                }
            }
        }

        y
    }

    /// Converts radians to degrees.
    ///
    /// # Examples
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
    #[inline]
    pub fn to_degrees(self) -> Self
    {
        self/(Self::FRAC_PI_2()/Self::from(90.0))
    }
    
    /// Converts degrees to radians.
    ///
    /// # Examples
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
    #[inline]
    pub fn to_radians(self) -> Self
    {
        self*(Self::FRAC_PI_2()/Self::from(90.0))
    }

    /// Returns the maximum of the two numbers.
    ///
    /// # Examples
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
    #[inline]
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
    /// # Examples
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
    #[inline]
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
    /// # Examples
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
                known as `fdimf` in C)."
    )]
    #[inline]
    pub fn abs_sub(self, other: Self) -> Self
    {
        if other >= self
        {
            return Self::zero()
        }
        self - other
    }

    /// Take the cubic root of a number.
    ///
    /// # Examples
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
    #[inline]
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
            let xabs = self.abs();
            if !SIGN_BIT && xabs < Self::one()
            {
                let xabs_log = xabs.recip().logb();
    
                let n_xabs_log = xabs_log/Self::from(3.0);
                
                n_xabs_log.expb().recip().copysign(self)
            }
            else
            {
                let xabs_log = xabs.logb();

                let n_xabs_log = xabs_log/Self::from(3.0);
                
                n_xabs_log.expb().copysign(self)
            }
        };

        const NEWTON: usize = NEWTON_RT;
        let third = Self::from(3.0).recip();
        let two = Self::from(2.0);
        let mut y = y;
        for _ in 0..NEWTON
        {
            let y_ = third*(self/(y*y) + two*y);
            if !y_.is_finite()
            {
                break
            }
            y = y_;
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
    #[inline]
    pub fn hypot(self, other: Self) -> Self
    {
        let bias = Self::exp_bias();
        let hi = Self::from(bias.to_f64().unwrap()*0.7).expb();
        let lo = Self::from(-bias.to_f64().unwrap()*0.7).expb();
    
        let mut x = self.abs();
        let mut y = other.abs();
    
        /* arrange |x| >= |y| */
        if x < y
        {
            core::mem::swap(&mut x, &mut y)
        }
    
        /* special cases */
        let ex = x.exp_bits();
        let ey = y.exp_bits();
        /* note: hypot(inf,nan) == inf */
        if y.is_infinite() || x.is_zero()
        {
            return y;
        }
        if x.is_infinite() || y.is_zero()
        {
            return x;
        }
        /* note: hypot(x,y) ~= x + y*y/x/2 with inexact for small y/x */
        if match <usize as NumCast>::from(ex - ey)
        {
            Some(de) => de > util::count_digits_in_base(FRAC_SIZE + 1, EXP_BASE),
            None => true
        }
        {
            return x + y;
        }
    
        /* precise sqrt argument in nearest rounding mode without overflow */
        /* xh*xh must not overflow and xl*xl must not underflow in sq */
        let bias_half = bias/U::from(2).unwrap();
        let mut z = Self::one();
        if ex > bias + bias_half
        {
            z = hi;
            x *= lo;
            y *= lo;
        }
        else if ey < bias - bias_half
        {
            z = lo;
            x *= hi;
            y *= hi;
        }
        z * (x*x + y*y).sqrt()
    }

    fn sin_extra_sign(self) -> (Self, bool)
    {
        if self.is_nan()
        {
            return (self, false)
        }
        if self.is_infinite()
        {
            return (Self::snan(), false)
        }
        if self.abs() < Self::from(0.000244140625)
        {
            return (self, false)
        }

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
        let mut p = unsafe {
            #[allow(static_mut_refs)]
            *P.get_or_insert_with(|| util::chebychev_approximation(C))
        };

        let one = Self::one();
        let mut w = self*Self::FRAC_2_PI();
        let mut i = 0;
        let four = Self::from(4.0);
        w %= four;
        while i < 4
        {
            w += four - one;
            if i % 2 == 0 && w < four
            {
                w = four - w;
            }
            w %= four;
            i += 1;
        }
        let two = Self::from(2.0);
        let mut s = false;
        let w = if w > one
        {
            if !SIGN_BIT && two < w
            {
                s = !s;
                w - two
            }
            else
            {
                two - w
            }
        }
        else if w < -one
        {
            -two - w
        }
        else
        {
            w
        };
        
        if w.abs() < Self::from(1.5542474911317903883680055016847e-4)
        {
            return (w*Self::FRAC_PI_2(), false)
        }

        let ww2 = two*w*w;
        let y = if !SIGN_BIT && ww2 < one
        {
            let z_neg = one - ww2;
            
            let mut y = w*util::polynomial(&p, z_neg, SIGN_BIT, true);
            if !SIGN_BIT && y.is_nan()
            {
                for p in p.iter_mut()
                {
                    *p = -*p
                }
                s = !s;
                y = w*util::polynomial(&p, z_neg, SIGN_BIT, true)
            }

            y
        }
        else
        {
            let z = ww2 - one;

            let mut y = w*util::polynomial(&p, z, SIGN_BIT, false);
            if !SIGN_BIT && y.is_nan()
            {
                for p in p.iter_mut()
                {
                    *p = -*p
                }
                s = !s;
                y = w*util::polynomial(&p, z, SIGN_BIT, false)
            }

            y
        };

        (y, s)
    }

    /// Computes the sine of a number (in radians).
    ///
    /// This implementation is based on Harvey M. Wagner's [Polynomial approximations to elementary functions](https://www.ams.org/journals/mcom/1954-08-047/S0025-5718-1954-0063487-2/S0025-5718-1954-0063487-2.pdf).
    ///
    /// # Examples
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
    #[inline]
    pub fn sin(self) -> Self
    {
        let (mut y, s) = self.sin_extra_sign();
        if s
        {
            y = -y
        }
        y
    }

    fn cos_extra_sign(self) -> (Self, bool)
    {
        if self.is_nan()
        {
            return (self, false)
        }
        if self.is_infinite()
        {
            return (Self::snan(), false)
        }
        
        const N: usize = 6;
        const C: [f64; N] = [
            0.472001216,
            -0.499403258,
            0.027992080,
            -0.000596695,
            0.000006704,
            -0.000000047
        ];

        static mut P: Option<[f64; N]> = None;
        let p = unsafe {
            #[allow(static_mut_refs)]
            *P.get_or_insert_with(|| util::chebychev_approximation(C))
        };

        let one = Self::one();
        let mut w = self*Self::FRAC_2_PI();
        let mut s = false;
        let mut i = 0;
        let four = Self::from(4u8);
        w %= four;
        while i < 4
        {
            w += four - one;
            if i % 2 == 1 && w < four
            {
                w = four - w;
            }
            w %= four;
            i += 1;
        }
        if w.abs() > one
        {
            s = !s
        }
        let two = Self::from(2u8);
        w = w.abs();
        let mut b = false;
        let w = if w > one
        {
            if !SIGN_BIT && two < w
            {
                b = !b;
                s = !s;
                w - two
            }
            else
            {
                two - w
            }
        }
        else
        {
            w
        };

        let p_neg = p.map(Neg::neg);

        let ww2 = two*w*w;
        let y = if !SIGN_BIT && ww2 < one
        {
            let z_neg = one - ww2;

            let mut y = util::polynomial(&if s {p_neg} else {p}, z_neg, SIGN_BIT, true);
            if !SIGN_BIT && y.is_nan()
            {
                b = !b;
                y = util::polynomial(&if s {p} else {p_neg}, z_neg, SIGN_BIT, true);
            }
            
            y
        }
        else
        {
            let z = ww2 - one;

            let mut y = util::polynomial(&if !SIGN_BIT && s {p_neg} else {p}, z, SIGN_BIT, false);
            if !SIGN_BIT && y.is_nan()
            {
                b = !b;
                y = util::polynomial(&if s {p} else {p_neg}, z, SIGN_BIT, false);
            }
            if SIGN_BIT && s
            {
                y = -y
            }

            y
        };

        (y, b)
    }

    /// Computes the cosine of a number (in radians).
    ///
    /// This implementation is based on Harvey M. Wagner's [Polynomial approximations to elementary functions](https://www.ams.org/journals/mcom/1954-08-047/S0025-5718-1954-0063487-2/S0025-5718-1954-0063487-2.pdf).
    ///
    /// # Examples
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
    /// assert!(abs_difference < FpDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn cos(self) -> Self
    {
        let (mut y, s) = self.cos_extra_sign();
        if s
        {
            y = -y
        }
        y
    }

    /// Computes the tangent of a number (in radians).
    ///
    /// # Examples
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
    /// assert!(abs_difference < FpDouble::from(1e-9));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn tan(mut self) -> Self
    {
        self %= Self::PI();
        while self > Self::FRAC_PI_2()
        {
            self -= Self::PI()
        }
        while self < -Self::FRAC_PI_2()
        {
            self += Self::PI()
        }

        let (sin, cos) = self.sin_cos();

        if sin.is_nan()
        {
            return sin
        }
        if cos.is_nan()
        {
            return cos
        }

        let mut y = sin/cos;

        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_TRIG;

            for _ in 0..NEWTON
            {
                let x = y.atan();
                let (dx, s) = if SIGN_BIT || x >= self
                {
                    (x - self, false)
                }
                else
                {
                    (self - x, true)
                };
                let dy = dx*(y*y + Self::one());
                if !dy.is_finite() || (!SIGN_BIT && !s && dy > y)
                {
                    break
                }
                if s
                {
                    y += dy
                }
                else
                {
                    y -= dy
                }
            }
        }

        if y.is_nan()
        {
            return Self::infinity().copysign(self)
        }

        y
    }

    /// Computes the arcsine of a number. Return value is in radians in
    /// the range [-pi/2, pi/2] or NaN if the number is outside the range
    /// [-1, 1].
    ///
    /// This implementation is based on Harvey M. Wagner's [Polynomial approximations to elementary functions](https://www.ams.org/journals/mcom/1954-08-047/S0025-5718-1954-0063487-2/S0025-5718-1954-0063487-2.pdf).
    ///
    /// # Examples
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
    /// assert!(abs_difference < FpDouble::from(1e-4));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn asin(self) -> Self
    {
        if self.is_nan()
        {
            return self
        }
        let xabs = self.abs();
        if xabs.is_one()
        {
            return Self::FRAC_PI_2().copysign(self)
        }
        if xabs > Self::one()
        {
            return Self::snan()
        }
        if xabs < Self::from(0.000000007450580596923828125)
        {
            return self
        }
        
        const N: usize = 10;
        const C: [f64; N] = [
            1.051231959,
            0.054946487,
            0.004080631,
            0.000407890,
            0.000046985,
            0.000005881,
            0.000000777,
            0.000000107,
            0.000000015,
            0.000000002
        ];

        static mut P: Option<[f64; N]> = None;
        let p = unsafe {
            #[allow(static_mut_refs)]
            P.get_or_insert_with(|| util::chebychev_approximation(C))
        };

        let one = Self::one();
        let w = if xabs <= Self::FRAC_1_SQRT_2()
        {
            self.abs()
        }
        else
        {
            (one - self*self).sqrt()
        };
        
        let ww2 = Self::from_uint(2u8)*w*w;
        
        let mut s = false;
        let mut y = if !SIGN_BIT && ww2 < one
        {
            let z_neg = one - ww2;

            let mut y = util::polynomial(p, z_neg, SIGN_BIT, true);
            if y.is_nan()
            {
                s = !s;
                y = util::polynomial(&p.map(Neg::neg), z_neg, SIGN_BIT, true)
            }
            y
        }
        else
        {
            let z = ww2 - one;

            util::polynomial(p, z, SIGN_BIT, false)
        }*w;
        if xabs > Self::FRAC_1_SQRT_2()
        {
            if s
            {
                y = Self::FRAC_PI_2() + y
            }
            else
            {
                y = Self::FRAC_PI_2() - y
            }
        }
        else if s
        {
            y = -y
        }
        if self.is_sign_negative()
        {
            y = -y
        }

        /*const PIO2_HI: f64 = 1.57079637050628662109375;
        const PIO2_LO: f64 = -4.37113900018624283e-8;
        const PIO4_HI: f64 = 0.785398185253143310546875;

        const P0: f64 = 1.666675248e-1;
        const P1: f64 = 7.495297643e-2;
        const P2: f64 = 4.547037598e-2;
        const P3: f64 = 2.417951451e-2;
        const P4: f64 = 4.216630880e-2;

        let pio2_hi = Self::from(PIO2_HI);
        let pio2_lo = Self::from(PIO2_LO);
        let pio4_hi = Self::from(PIO4_HI);

        let mut y = if xabs.is_one()
        {
            self*pio2_hi + self*pio2_lo
        }
        else if xabs > Self::one()
        {
            return (self - self)/(self - self)
        }
        else if xabs < Self::from(0.5)
        {
            if xabs < Self::from(0.000000007450580596923828125)
            {
                self
            }
            else
            {
                let t = self*self;
                let w = t*(Self::from(P0) + t*(Self::from(P1) + t*(Self::from(P2) + t*(Self::from(P3) + t*Self::from(P4)))));
                self + self*w
            }
        }
        else
        {
            let mut w = Self::one() - xabs;
            let mut t = w*Self::from(0.5);
            let mut p = t*(Self::from(P0) + t*(Self::from(P1) + t*(Self::from(P2) + t*(Self::from(P3) + t*Self::from(P4)))));
            let s = t.sqrt();

            let two = Self::from_uint(2u8);

            if xabs >= Self::from(0.975)
            {
                t = pio2_hi - (two*(s + s*p) - pio2_lo);
            }
            else
            {
                let e_s = s.exp_bits();
                let mut f_s = s.frac_bits();
                if !Self::IS_INT_IMPLICIT
                {
                    f_s = f_s + (s.int_bits() << Self::INT_POS)
                }
                f_s = f_s & ((U::max_value() >> FRAC_SIZE/2) << FRAC_SIZE/2);
                w = Self::from_bits((e_s << Self::EXP_POS) + (f_s << Self::FRAC_POS));

	            let c = (t - w*w)/(s + w);
                let r = p;
                p = two*s*r - (pio2_lo - two*c);
                let q  = pio4_hi - two*w;
                t = pio4_hi - (p - q);
            }

            t.copysign(self)
        };*/

        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_TRIG;

            for _ in 0..NEWTON
            {
                let (sin, sin_s) = y.sin_extra_sign();
                let (cos, cos_s) = y.cos_extra_sign();
                let (dsx, dsx_s) = util::add_extra_sign((sin, sin_s), (self, true), SIGN_BIT);
                let s = dsx_s^cos_s;
                let dy = dsx/cos;
                if !dy.is_finite() || (!SIGN_BIT && !s && dy > y)
                {
                    break
                }
                if s
                {
                    y += dy
                }
                else
                {
                    y -= dy
                }
            }
        }

        y
    }

    /// Computes the arccosine of a number. Return value is in radians in
    /// the range [0, pi] or NaN if the number is outside the range
    /// [-1, 1].
    ///
    /// This implementation is based on Harvey M. Wagner's [Polynomial approximations to elementary functions](https://www.ams.org/journals/mcom/1954-08-047/S0025-5718-1954-0063487-2/S0025-5718-1954-0063487-2.pdf).
    ///
    /// # Examples
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
        let one = Self::one();
        if self == one
        {
            return Self::zero()
        }
        if self == -one
        {
            return Self::PI()
        }
        let xabs = self.abs();
        if xabs > one
        {
            return Self::snan()
        }
        if xabs < Self::from(0.00000001490116119384765625)
        {
            return Self::FRAC_PI_2() - self
        }
        
        const N: usize = 10;
        const C: [f64; N] = [
            1.051231959,
            0.054946487,
            0.004080631,
            0.000407890,
            0.000046985,
            0.000005881,
            0.000000777,
            0.000000107,
            0.000000015,
            0.000000002
        ];

        static mut P: Option<[f64; N]> = None;
        let p = unsafe {
            #[allow(static_mut_refs)]
            P.get_or_insert_with(|| util::chebychev_approximation(C))
        }.map(Self::from);

        let w = if xabs <= Self::FRAC_1_SQRT_2()
        {
            self.abs()
        }
        else
        {
            (one - self*self).sqrt()
        };
        
        let ww2 = Self::from_uint(2u8)*w*w;
        let z = ww2 - Self::one();

        let mut y = util::polynomial(&p, z, true, false)*w;
        if xabs <= Self::FRAC_1_SQRT_2()
        {
            y = Self::FRAC_PI_2() - y
        }
        if self.is_sign_negative()
        {
            y = Self::PI() - y
        }
        
        if !y.is_finite()
        {
            let xx = self*self;
            y = if xx > Self::from(0.5)
            {
                ((Self::one() - xx).sqrt()/self).atan()
            }
            else
            {
                Self::FRAC_PI_2() - (self/(Self::one() - xx).sqrt()).atan()
            };
        }

        /*const PIO2_HI: f64 = 1.5707962513e+00;
        const PIO2_LO: f64 = 7.5497894159e-08;

        const PS0: f64 = 1.6666667163e-01;
        const PS1: f64 = -3.2556581497e-01;
        const PS2: f64 = 2.0121252537e-01;
        const PS3: f64 = -4.0055535734e-02;
        const PS4: f64 = 7.9153501429e-04;
        const PS5: f64 = 3.4793309169e-05;
        const QS1: f64 = -2.4033949375e+00;
        const QS2: f64 = 2.0209457874e+00;
        const QS3: f64 = -6.8828397989e-01;
        const QS4: f64 = 7.7038154006e-02;

        let mut y = if xabs < Self::from(0.5)
        {
            if xabs < Self::from(0.00000001490116119384765625)
            {
                Self::from(PIO2_HI) + Self::from(PIO2_LO)
            }
            else
            {
                let z = self*self;
                let p = z*(Self::from(PS0) + z*(Self::from(PS1) + z*(Self::from(PS2) + z*(Self::from(PS3) + z*(Self::from(PS4) + z*Self::from(PS5))))));
                let q = Self::one() + z*(Self::from(QS1) + z*(Self::from(QS2) + z*(Self::from(QS3) + z*Self::from(QS4))));
                let r = p/q;
                Self::from(PIO2_HI) - (self - (Self::from(PIO2_LO) - self*r))
            }
        }
        else if self.is_sign_negative()
        {
            let z = (Self::one() + self)*Self::from(0.5);
            let p = z*(Self::from(PS0) + z*(Self::from(PS1) + z*(Self::from(PS2) + z*(Self::from(PS3) + z*(Self::from(PS4) + z*Self::from(PS5))))));
            let q = Self::one() + z*(Self::from(QS1) + z*(Self::from(QS2) + z*(Self::from(QS3) + z*Self::from(QS4))));
            let r = p/q;
            let s = z.sqrt();
            let w = r*s - Self::from(PIO2_LO);
            Self::PI() - Self::from_uint(2u8)*(s + w)
        }
        else
        {
            let z = (Self::one() - self)*Self::from(0.5);
            let s = z.sqrt();

            let e_s = s.exp_bits();
            let mut f_s = s.frac_bits();
            if !Self::IS_INT_IMPLICIT
            {
                f_s = f_s + (s.int_bits() << Self::INT_POS)
            }
            f_s = f_s & ((U::max_value() >> FRAC_SIZE/2) << FRAC_SIZE/2);
            let df = Self::from_bits((e_s << Self::EXP_POS) + (f_s << Self::FRAC_POS));
            let c = (z - df*df)/(s + df);

            let p = z*(Self::from(PS0) + z*(Self::from(PS1) + z*(Self::from(PS2) + z*(Self::from(PS3) + z*(Self::from(PS4) + z*Self::from(PS5))))));
            let q = Self::one() + z*(Self::from(QS1) + z*(Self::from(QS2) + z*(Self::from(QS3) + z*Self::from(QS4))));
            let r = p/q;
            let w = r*s + c;
            Self::from_uint(2u8)*(df + w)
        };*/

        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_TRIG;

            for _ in 0..NEWTON
            {
                let (cos, cos_s) = y.cos_extra_sign();
                let (sin, sin_s) = y.sin_extra_sign();
                let (dcx, dcx_s) = util::add_extra_sign((cos, cos_s), (self, true), SIGN_BIT);
                let s = dcx_s^sin_s;
                let dy = dcx/sin;
                if !dy.is_finite() || (!SIGN_BIT && s && dy > y)
                {
                    break
                }
                if s
                {
                    y -= dy
                }
                else
                {
                    y += dy
                }
            }
        }

        y
    }

    /// Computes the arctangent of a number. Return value is in radians in the
    /// range [-pi/2, pi/2];
    ///
    /// This implementation is based on Harvey M. Wagner's [Polynomial approximations to elementary functions](https://www.ams.org/journals/mcom/1954-08-047/S0025-5718-1954-0063487-2/S0025-5718-1954-0063487-2.pdf).
    ///
    /// # Examples
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
    /// assert!(abs_difference < FpDouble::from(1e-5));
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
        
        const N: usize = 11;
        const C: [f64; N] = [
            0.881373587,
            -0.105892925,
            0.011135843,
            -0.001381195,
            0.000185743,
            -0.000026215,
            0.000003821,
            -0.00000057,
            0.000000086,
            -0.000000013,
            0.000000002
        ];

        static mut P: Option<[f64; N]> = None;
        let p = unsafe {
            #[allow(static_mut_refs)]
            P.get_or_insert_with(|| util::chebychev_approximation(C))
        };

        let xabs = self.abs();
        let inv = xabs > Self::one();
        let w = if !inv
        {
            self
        }
        else
        {
            self.recip()
        };
        
        let one = Self::one();
        let ww2 = Self::from_uint(2u8)*w*w;
        let mut s = false;
        let mut y = if !SIGN_BIT && ww2 < one
        {
            let z = one - ww2;

            let mut y = util::polynomial(p, z, SIGN_BIT, true);
            if y.is_nan()
            {
                s = !s;
                y = util::polynomial(&p.map(Neg::neg), z, SIGN_BIT, true);
            }
            y
        }
        else
        {
            let z = ww2 - one;

            let mut y = util::polynomial(p, z, SIGN_BIT, false);
            if !SIGN_BIT && y.is_nan()
            {
                s = !s;
                y = util::polynomial(&p.map(Neg::neg), z, SIGN_BIT, false);
            }
            y
        }*w;
        if inv
        {
            if s
            {
                y = Self::FRAC_PI_2() + y
            }
            else if !SIGN_BIT && y > Self::FRAC_PI_2()
            {
                y %= Self::PI();
                y = Self::PI() + Self::FRAC_PI_2() - y
            }
            else
            {
                y = Self::FRAC_PI_2() - y
            }
        }
        else if s
        {
            y %= Self::PI();
            y = Self::PI() - y
        }

        /*if y.is_nan()
        {
            const TAYLOR: usize = 8;
            y = if self.abs() < Self::one()
            {
                let mut s = false;
                let mut z = self;
                let mut y = z;

                for k in 1..TAYLOR
                {
                    z *= self*self;
                    s = !s;
                    let dy = z/Self::from_uint(1 + 2*k);
                    if !s
                    {
                        y += dy
                    }
                    else
                    {
                        y -= dy
                    }
                }

                y
            }
            else
            {
                let mut s = false;
                let mut z = Self::one()/self;
                let mut y = Self::FRAC_PI_2() - z;

                for k in 1..TAYLOR
                {
                    z /= self*self;
                    s = !s;
                    let dy = z/Self::from_uint(1 + 2*k);
                    if !s
                    {
                        y -= dy
                    }
                    else
                    {
                        y += dy
                    }
                }

                y
            }
        }*/

        y %= Self::PI();
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
                let (sin, sin_s) = y.sin_extra_sign();
                let (cos, cos_s) = y.cos_extra_sign();
                let xcos = self*cos;
                let (dsxc, dsxc_s) = util::add_extra_sign((sin, sin_s), (xcos, !cos_s), SIGN_BIT);
                let dy = dsxc*cos;
                let s = dsxc_s^cos_s;
                if !dy.is_finite() || (!SIGN_BIT && !s && dy > y)
                {
                    break
                }
                if s
                {
                    y += dy
                }
                else
                {
                    y -= dy
                }
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
    /// # Examples
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
    /// assert!(abs_difference_1 < FpDouble::from(1e-5));
    /// assert!(abs_difference_2 < FpDouble::from(1e-5));
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
    /// # Examples
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
        (self.sin(), self.cos())
    }

    /// Returns `e^(self) - 1`.
    ///
    /// # Examples
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
        if self > Self::one()
        {
            return self.exp() - Self::one()
        }
        if self.is_zero()
        {
            return Self::zero()
        }

        let mut y = self.exp_nonewton() - Self::one();

        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_EXP;

            for _ in 0..NEWTON
            {
                let yp1 = y + Self::one();
                let x = y.ln_1p();
                let (dy, s) = if SIGN_BIT || x >= self
                {
                    (yp1*(x - self), false)
                }
                else
                {
                    (yp1*(self - x), true)
                };
                if !dy.is_finite() || (!SIGN_BIT && !s && dy > y)
                {
                    break
                }
                if s
                {
                    y += dy
                }
                else
                {
                    y -= dy
                }
            }
        }

        y
    }

    /// Returns `ln(1+n)` (natural logarithm) more accurately than if
    /// the operations were performed separately.
    ///
    /// # Examples
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

            let one = Self::one();
            for _ in 0..NEWTON
            {
                let xp1dx = xp1/y.exp();
                let (dy, s) = if SIGN_BIT || one >= xp1dx
                {
                    (one - xp1dx, false)
                }
                else
                {
                    (xp1dx - one, true)
                };
                if !dy.is_finite() || (!SIGN_BIT && !s && dy > y)
                {
                    break
                }
                if s
                {
                    y += dy
                }
                else
                {
                    y -= dy
                }
            }
        }

        y
    }

    /// Hyperbolic sine function.
    ///
    /// # Examples
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
    /// assert!(abs_difference < FpDouble::from(1e-3));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn sinh(self) -> Self
    {
        if !self.is_finite()
        {
            return self
        }

        let emx = (-self.abs()).exp();
    
        let mut y = ((Self::one() - emx*emx)/emx*Self::from(0.5)).copysign(self);

        if y.is_nan()
        {
            let ex = (self.abs()).exp();
        
            y = ((ex*ex - Self::one())/ex*Self::from(0.5)).copysign(self);
        }

        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_TRIG;

            for _ in 0..NEWTON
            {
                let x = y.asinh();
                let (dy, s) = if SIGN_BIT || x >= self
                {
                    ((x - self)*(y*y + Self::one()).sqrt(), false)
                }
                else
                {
                    ((self - x)*(y*y + Self::one()).sqrt(), true)
                };
                if !dy.is_finite() || (!SIGN_BIT && !s && dy > y)
                {
                    break
                }
                if s
                {
                    y += dy
                }
                else
                {
                    y -= dy
                }
            }
        }

        if y.is_nan()
        {
            return Self::infinity().copysign(self)
        }

        y
    }

    /// Hyperbolic cosine function.
    ///
    /// # Examples
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
    /// assert!(abs_difference < FpDouble::from(1.0e-3));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn cosh(self) -> Self
    {
        let xabs = self.abs();
        if !self.is_finite()
        {
            return xabs
        }

        let emx = (-xabs).exp();
    
        let mut y = (Self::one() + emx*emx)/emx*Self::from(0.5);

        if !y.is_finite()
        {
            let ex = xabs.exp();
        
            y = (ex*ex + Self::one())/ex*Self::from(0.5);
        }

        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_TRIG;

            for _ in 0..NEWTON
            {
                let x = y.acosh();
                let (dy, s) = if SIGN_BIT || x >= xabs
                {
                    ((x - xabs)*(y*y - Self::one()).sqrt(), false)
                }
                else
                {
                    ((xabs - x)*(y*y - Self::one()).sqrt(), true)
                };
                if !dy.is_finite() || (!SIGN_BIT && !s && dy > y)
                {
                    break
                }
                if s
                {
                    y += dy
                }
                else
                {
                    y -= dy
                }
            }
        }

        y
    }
    
    /// Hyperbolic tangent function.
    ///
    /// # Examples
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
        let one = Self::one();
        let xabs = self.abs();
        if xabs < Self::from(2.7755575615628913510590791702271e-17)
        {
            return self*(one + self)
        }

        let mut y = if SIGN_BIT
        {
            let ex = (-xabs).exp();
            let ex2 = ex*ex;
            let ex2p1 = one + ex2;
            let ex2m1 = one - ex2;
            
            (ex2m1/ex2p1).copysign(self)
        }
        else
        {
            let ex = xabs.exp().recip();
            let ex2 = ex*ex;
            let ex2p1 = one + ex2;
            let ex2m1 = one - ex2;
            
            ex2m1/ex2p1
        };

        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_TRIG;

            for _ in 0..NEWTON
            {
                let x = y.atanh();
                if !x.is_finite()
                {
                    break
                }
                let (dy, s) = if SIGN_BIT || x >= self
                {
                    ((x - self)*(one - y)*(one + y), false)
                }
                else
                {
                    ((self - x)*(one - y)*(one + y), true)
                };
                if !dy.is_finite() || (!SIGN_BIT && !s && dy > y)
                {
                    break
                }
                if s
                {
                    y += dy
                }
                else
                {
                    y -= dy
                }
            }
        }

        if y.abs() >= Self::one()
        {
            return Self::one().copysign(y)
        }
        if y.is_nan()
        {
            return Self::one().copysign(self)
        }

        y
    }

    /// Inverse hyperbolic sine function.
    ///
    /// This implementation is based on the glibc implementation of asinhf.
    ///
    /// # Examples
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
        //(self + (self*self + Self::one()).sqrt()).ln()

        let xabs = self.abs();
        if xabs < Self::from(0.00006103515625)
        {
            return self
        }
        let w = if xabs > Self::from(16384.0)
        {
            if !xabs.is_finite()
            {
                return self + self
            }
            xabs.ln() + Self::LN_2()
        }
        else
        {
            let one = Self::one();
            let two = Self::from_uint(2u8);
            let t = xabs*xabs;
            if xabs > two
            {
                (two*xabs + ((t + one).sqrt() + xabs).recip()).ln()
            }
            else
            {
                (xabs + t/(one + (one + t).sqrt())).ln_1p()
            }
        };
        w.copysign(self)
    }

    /// Inverse hyperbolic cosine function.
    /// 
    /// This implementation is based on the glibc implementation of acoshf.
    ///
    /// # Examples
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
    /// assert!(abs_difference < FpDouble::from(1.0e-3));
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
            return Self::snan()
        }
        if self.is_infinite()
        {
            return self
        }
        //(self + (self*self - Self::one()).sqrt()).ln()

        if self > Self::from(268435456.0)
        {
            self.ln() + Self::LN_2()
        }
        else if self.is_one()
        {
            Self::zero()
        }
        else
        {
            let two = Self::from_uint(2u8);
            let one = Self::one();
            if self > two
            {
                let t = self*self;
                (two*self - (self + (t - one).sqrt()).recip()).ln()
            }
            else
            {
                let t = self - one;
                (t + (two*t + t*t).sqrt()).ln_1p()
            }
        }
    }

    /// Inverse hyperbolic tangent function.
    ///
    /// This implementation is based on the glibc implementation of atanhf.
    ///
    /// # Examples
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
        let xabs = self.abs();
        if xabs > Self::one()
        {
            return Self::snan()
        }
        //Self::from(0.5)*((Self::one() + self.abs())/(Self::one() - self.abs())).ln().copysign(self)
        
        let mut t = xabs + xabs;
        let one = Self::one();
        let half = Self::from(0.5);
        if xabs < half
        {
            if xabs < Self::from(0.0000000037252902984619140625)
            {
                return self
            }

            t = half*(t + t*xabs/(one - xabs)).ln_1p()
        }
        else
        {
            t = half*(t/(one - xabs)).ln_1p()
        }

        t.copysign(self)
    }
    
    fn ln_gamma_lanczos(self) -> Self
    {
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

        let mut sum = Self::from(LANCZOS_CHEB_7[0]);
        for (k, c) in LANCZOS_CHEB_7.iter()
            .skip(1)
            .enumerate()
        {
            if SIGN_BIT || c.is_sign_positive()
            {
                sum += Self::from(*c)/(self + Self::from_uint(k))
            }
        }
        if !SIGN_BIT
        {
            for (k, c) in LANCZOS_CHEB_7.iter()
                .skip(1)
                .enumerate()
            {
                if c.is_sign_negative()
                {
                    sum -= Self::from(-c)/(self + Self::from_uint(k))
                }
            }
        }

        let term1 = (self - Self::from(0.5))
            *((self + Self::from(6.5))/Self::E()).ln();
        let term2 = Self::from(LOGROOT2PI) + sum.ln();

        let seven = Self::from(7u8);
        if term2 < seven
        {
            return term1 + term2 - seven
        }
        term1 + (term2 - seven)
    }

    /// Natural logarithm of the absolute value of the gamma function
    ///
    /// The integer part of the tuple indicates the sign of the gamma function.
    ///
    /// This implementation is based on [the libstdc++ implementation of the gamma function](https://gcc.gnu.org/onlinedocs/libstdc++/libstdc++-html-USERS-4.4/a01203.html).
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
        if !self.is_finite()
        {
            return (self.abs(), if self.is_nan() || self.is_sign_negative() {0} else {1})
        }
        if self.is_zero()
        {
            return (Self::infinity(), if self.is_sign_negative() {-1} else {1})
        }

        if self >= Self::from(0.5)
        {
            return (self.ln_gamma_lanczos(), 1)
        }

        let sin_fact = (Self::PI()*self).sin();

        if sin_fact.is_zero()
        {
            return (Self::infinity(), 0)
        }

        let sin_fact_abs = sin_fact.abs();

        (
            if !SIGN_BIT && sin_fact_abs < Self::one()
            {
                Self::PI().ln() + sin_fact_abs.recip().ln() - (Self::one() - self).ln_gamma_lanczos()
            }
            else
            {
                Self::PI().ln() - sin_fact_abs.ln() - (Self::one() - self).ln_gamma_lanczos()
            },
            if sin_fact.is_sign_negative() {-1} else {1}
        )
    }

    /// Gamma function.
    /// 
    /// This implementation is based on [the libstdc++ implementation of the gamma function](https://gcc.gnu.org/onlinedocs/libstdc++/libstdc++-html-USERS-4.4/a01203.html).
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
        if !SIGN_BIT
        {
            return self
        }

        let mask = (U::max_value() >> (util::bitsize_of::<U>() - Self::SIGN_SIZE)) << Self::SIGN_POS;
        
        Self::from_bits((self.to_bits() & (!mask)) | (sign.to_bits() & mask))
    }

    fn with_sign(self, neg: bool) -> Self
    {
        if !SIGN_BIT
        {
            if !neg || self.is_zero()
            {
                return self
            }
            return Self::qnan()
        }
        
        let mask = (U::max_value() >> (util::bitsize_of::<U>() - Self::SIGN_SIZE)) << Self::SIGN_POS;

        let mut bits = self.to_bits() & (!mask);
        if neg
        {
            bits = bits | mask;
        }
        
        Self::from_bits(bits)
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
            return Self::from_bits(Self::shift_exp(bias) + Self::shift_int(U::one()))
        }
        Self::from_bits(Self::shift_exp(bias))
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
    ///
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
    pub fn total_cmp(self, other: Self) -> Ordering
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
                3u8
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
            f0 = f0 + Self::shift_int(self.int_bits());
            f1 = f1 + Self::shift_int(other.int_bits());
            
            Self::normalize_mantissa_down(&mut e0, &mut f0, Some(e1));
            Self::normalize_mantissa_down(&mut e1, &mut f1, Some(e0));
            Self::normalize_mantissa_up(&mut e0, &mut f0, Some(e1));
            Self::normalize_mantissa_up(&mut e1, &mut f1, Some(e0));
        }
        
        if e0 != e1
        {
            return if s {e1.cmp(&e0)} else {e0.cmp(&e1)}
        }

        if s {f1.cmp(&f0)} else {f0.cmp(&f1)}
    }

    /// Returns `self*EXP_BASE`.
    ///
    /// This is generally faster than using regular multiplication.
    ///
    /// # Examples
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
    /// let abs_difference_f = (f.mulb() - FpDouble::from(4.0)).abs();
    ///
    /// // 2*10 - 20 == 0
    /// let abs_difference_d = (d.mulb() - DecDouble::from(20.0)).abs();
    ///
    /// assert!(abs_difference_f < FpDouble::from(1e-10));
    /// assert!(abs_difference_d < DecDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn mulb(self) -> Self
    {
        if !self.is_finite()
        {
            return self
        }
        let s = self.is_sign_negative();
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

        Self::normalize_mantissa(&mut e, &mut f, None);
        Self::from_sign_exp_mantissa(s, e, f)
    }

    /// Returns `self/EXP_BASE`.
    ///
    /// This is generally faster than using regular division.
    ///
    /// # Examples
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
    /// let abs_difference_f = (f.divb() - FpDouble::from(1.0)).abs();
    ///
    /// // 2/10 - 0.2 == 0
    /// let abs_difference_d = (d.divb() - DecDouble::from(0.2)).abs();
    ///
    /// assert!(abs_difference_f < FpDouble::from(1e-10));
    /// assert!(abs_difference_d < DecDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn divb(self) -> Self
    {
        if !self.is_finite()
        {
            return self
        }
        let s = self.is_sign_negative();
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
        
        Self::normalize_mantissa(&mut e, &mut f, None);
        Self::from_sign_exp_mantissa(s, e, f)
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

    fn erfc1(self) -> Self
    {
        const ERX: f64 = 8.45062911510467529297e-01; /* 0x3FEB0AC1, 0x60000000 */
        
        /*
        * Coefficients for approximation to  erf  in [0.84375,1.25]
        */
        const PA: [f64; 7] = [
            -2.36211856075265944077e-03, /* 0xBF6359B8, 0xBEF77538 */
            4.14856118683748331666e-01, /* 0x3FDA8D00, 0xAD92B34D */
            -3.72207876035701323847e-01, /* 0xBFD7D240, 0xFBB8C3F1 */
            3.18346619901161753674e-01, /* 0x3FD45FCA, 0x805120E4 */
            -1.10894694282396677476e-01, /* 0xBFBC6398, 0x3D3E28EC */
            3.54783043256182359371e-02, /* 0x3FA22A36, 0x599795EB */
            -2.16637559486879084300e-03 /* 0xBF61BF38, 0x0A96073F */
        ];
        const QA: [f64; 6] = [
            1.06420880400844228286e-01, /* 0x3FBB3E66, 0x18EEE323 */
            5.40397917702171048937e-01, /* 0x3FE14AF0, 0x92EB6F33 */
            7.18286544141962662868e-02, /* 0x3FB2635C, 0xD99FE9A7 */
            1.26171219808761642112e-01, /* 0x3FC02660, 0xE763351F */
            1.36370839120290507362e-02, /* 0x3F8BEDC2, 0x6B51DD1C */
            1.19844998467991074170e-02 /* 0x3F888B54, 0x5735151D */
        ];
        let one = Self::one();
        let xabs = self.abs();
        let (s, b) = if xabs >= one
        {
            (xabs - one, false)
        }
        else
        {
            (one - xabs, true)
        };
        let mut r = false;
        let mut p = util::polynomial(&PA, s, SIGN_BIT, b);
        if !SIGN_BIT && p.is_nan()
        {
            p = util::polynomial(&PA.map(Neg::neg), s, SIGN_BIT, b);
            r = !r;
        }
        let q = if b
        {
            let mut q = s*util::polynomial(&QA, s, SIGN_BIT, b);
            if !SIGN_BIT && q.is_nan()
            {
                q = s*util::polynomial(&QA.map(Neg::neg), s, SIGN_BIT, b);
                one + q
            }
            else if SIGN_BIT || one >= q
            {
                one - q
            }
            else
            {
                r = !r;
                q - one
            }
        }
        else
        {
            let mut q = s*util::polynomial(&QA, s, SIGN_BIT, b);
            if SIGN_BIT || !q.is_nan()
            {
                one + q
            }
            else
            {
                q = s*util::polynomial(&QA.map(Neg::neg), s, SIGN_BIT, b);
                if one >= q
                {
                    one - q
                }
                else
                {
                    r = !r;
                    q - one
                }
            }
        };
    
        if r
        {
            one - Self::from(ERX) + p/q
        }
        else
        {
            one - Self::from(ERX) - p/q
        }
    }
    
    fn erfc2(mut self) -> Self
    {
        /*
         * Coefficients for approximation to  erfc in [1.25,1/0.35]
         */
        const RA: [f64; 8] = [
            -9.86494403484714822705e-03, /* 0xBF843412, 0x600D6435 */
            -6.93858572707181764372e-01, /* 0xBFE63416, 0xE4BA7360 */
            -1.05586262253232909814e+01, /* 0xC0251E04, 0x41B0E726 */
            -6.23753324503260060396e+01, /* 0xC04F300A, 0xE4CBA38D */
            -1.62396669462573470355e+02, /* 0xC0644CB1, 0x84282266 */
            -1.84605092906711035994e+02, /* 0xC067135C, 0xEBCCABB2 */
            -8.12874355063065934246e+01, /* 0xC0545265, 0x57E4D2F2 */
            -9.81432934416914548592e+00 /* 0xC023A0EF, 0xC69AC25C */
        ];
        const SA: [f64; 8] = [
            1.96512716674392571292e+01, /* 0x4033A6B9, 0xBD707687 */
            1.37657754143519042600e+02, /* 0x4061350C, 0x526AE721 */
            4.34565877475229228821e+02, /* 0x407B290D, 0xD58A1A71 */
            6.45387271733267880336e+02, /* 0x40842B19, 0x21EC2868 */
            4.29008140027567833386e+02, /* 0x407AD021, 0x57700314 */
            1.08635005541779435134e+02, /* 0x405B28A3, 0xEE48AE2C */
            6.57024977031928170135e+00, /* 0x401A47EF, 0x8E484A93 */
            -6.04244152148580987438e-02 /* 0xBFAEEFF2, 0xEE749A62 */
        ];
        /*
         * Coefficients for approximation to  erfc in [1/.35,28]
         */
        const RB: [f64; 7] = [
            -9.86494292470009928597e-03, /* 0xBF843412, 0x39E86F4A */
            -7.99283237680523006574e-01, /* 0xBFE993BA, 0x70C285DE */
            -1.77579549177547519889e+01, /* 0xC031C209, 0x555F995A */
            -1.60636384855821916062e+02, /* 0xC064145D, 0x43C5ED98 */
            -6.37566443368389627722e+02, /* 0xC083EC88, 0x1375F228 */
            -1.02509513161107724954e+03, /* 0xC0900461, 0x6A2E5992 */
            -4.83519191608651397019e+02 /* 0xC07E384E, 0x9BDC383F */
        ];
        const SB: [f64; 7] = [
            3.03380607434824582924e+01, /* 0x403E568B, 0x261D5190 */
            3.25792512996573918826e+02, /* 0x40745CAE, 0x221B9F0A */
            1.53672958608443695994e+03, /* 0x409802EB, 0x189D5118 */
            3.19985821950859553908e+03, /* 0x40A8FFB7, 0x688C246A */
            2.55305040643316442583e+03, /* 0x40A3F219, 0xCEDF3BE6 */
            4.74528541206955367215e+02, /* 0x407DA874, 0xE79FE763 */
            -2.24409524465858183362e+01 /* 0xC03670E2, 0x42712D62 */
        ];

        if self.abs() < Self::from(1.25)
        {
            /* |x| < 1.25 */
            return self.erfc1();
        }

        self = self.abs();
        let mut b = false;
        let s = (self*self).recip();
        let mut r;
        let mut big_s;
        let one = Self::one();
        if self < Self::from(1.0/0.35)
        {
            /* |x| < 1/.35 ~ 2.85714 */
            r = util::polynomial(&RA, s, SIGN_BIT, false);
            if !SIGN_BIT && r.is_nan()
            {
                r = util::polynomial(&RA.map(Neg::neg), s, SIGN_BIT, false);
                b = !b;
            }
            big_s = one + s*util::polynomial(&SA, s, SIGN_BIT, false);
            if !SIGN_BIT && big_s.is_nan()
            {
                big_s = s*util::polynomial(&SA.map(Neg::neg), s, SIGN_BIT, false) - one;
                b = !b;
            }
        }
        else
        {
            /* |x| > 1/.35 */
            r = util::polynomial(&RB, s, SIGN_BIT, false);
            if !SIGN_BIT && r.is_nan()
            {
                r = util::polynomial(&RB.map(Neg::neg), s, SIGN_BIT, false);
                b = !b;
            }
            big_s = one + s*util::polynomial(&SB, s, SIGN_BIT, false);
            if !SIGN_BIT && big_s.is_nan()
            {
                big_s = s*util::polynomial(&SB.map(Neg::neg), s, SIGN_BIT, false) - one;
                b = !b;
            }
        }
        let z = Self::from_bits((self.to_bits() >> (FRAC_SIZE/2)) << (FRAC_SIZE/2));

        let ofs = Self::from(0.5625);
        if !SIGN_BIT
        {
            match z >= self
            {
                true => if b
                {
                    let p = (self - z) * (z + self);
                    let q = r / big_s;
                    if p >= q
                    {
                        (z * z + ofs).exp().recip() * (p - q).exp() / self
                    }
                    else
                    {
                        (z * z + ofs).exp().recip() / (q - p).exp() / self
                    }
                }
                else
                {
                    (z * z + ofs).exp().recip() * ((z - self) * (z + self) + r / big_s).exp() / self
                },
                false => if b
                {
                    (z * z + ofs).exp().recip() / ((z - self) * (z + self) + r / big_s).exp() / self
                }
                else
                {
                    let p = r / big_s;
                    let q = (self - z) * (z + self);
                    if p >= q
                    {
                        (z * z + ofs).exp().recip() * (p - q).exp() / self
                    }
                    else
                    {
                        (z * z + ofs).exp().recip() / (q - p).exp() / self
                    }
                }
            }
        }
        else
        {
            (-z * z - ofs).exp() * ((z - self) * (z + self) + r / big_s).exp() / self
        }
    }
        
    /// Error function (f64)
    ///
    /// Calculates an approximation to the “error function”, which estimates
    /// the probability that an observation will fall within x standard
    /// deviations of the mean (assuming a normal distribution)
    pub fn erf(self) -> Self
    {
        /*
        * Coefficients for approximation to  erf on [0,0.84375]
        */
        const EFX8: f64 = 1.02703333676410069053e+00; /* 0x3FF06EBA, 0x8214DB69 */
        const PP: [f64; 5] = [
            1.28379167095512558561e-01, /* 0x3FC06EBA, 0x8214DB68 */
            -3.25042107247001499370e-01, /* 0xBFD4CD7D, 0x691CB913 */
            -2.84817495755985104766e-02, /* 0xBF9D2A51, 0xDBD7194F */
            -5.77027029648944159157e-03, /* 0xBF77A291, 0x236668E4 */
            -2.37630166566501626084e-05 /* 0xBEF8EAD6, 0x120016AC */
        ];
        const QQ: [f64; 5] = [
            3.97917223959155352819e-01, /* 0x3FD97779, 0xCDDADC09 */
            6.50222499887672944485e-02, /* 0x3FB0A54C, 0x5536CEBA */
            5.08130628187576562776e-03, /* 0x3F74D022, 0xC4D36B0F */
            1.32494738004321644526e-04, /* 0x3F215DC9, 0x221C1A10 */
            -3.96022827877536812320e-06 /* 0xBED09C43, 0x42A26120 */
        ];

        let sign = self.sign_bit();
        let xabs = self.abs();
        let one = Self::one();
        if self.is_nan()
        {
            return self
        }
        if self.is_infinite()
        {
            return Self::one().copysign(self)
        }
        let y;
        if xabs < Self::from(0.84375)
        {
            /* |x| < 0.84375 */
            if xabs < Self::from(0.0000000037252902984619140625)
            {
                /* |x| < 2**-28 */
                /* avoid underflow */
                return Self::from(0.125)*(Self::from(8u8)*self + Self::from(EFX8) * self);
            }
            let z = self*self;
            let mut b = false;
            let mut r = util::polynomial(&PP, z, SIGN_BIT, false);
            let mut s = one + z*util::polynomial(&QQ, z, SIGN_BIT, false);
            if !SIGN_BIT && r.is_nan()
            {
                r = util::polynomial(&PP.map(Neg::neg), z, SIGN_BIT, false);
                b = !b;
            }
            if !SIGN_BIT && s.is_nan()
            {
                s = z*util::polynomial(&QQ.map(Neg::neg), z, SIGN_BIT, false) - one;
                b = !b;
            }
            y = r / s;
            if b
            {
                return self - self*y;
            }
            return self + self*y;
        }
        if xabs < Self::from_uint(6u8)
        {
            /* 0.84375 <= |x| < 6 */
            y = one - self.erfc2();
        }
        else
        {
            let x1p_1022 = Self::from(f64::from_bits(0x0010000000000000));
            y = one - x1p_1022;
        }
    
        if !sign.is_zero()
        {
            -y
        }
        else
        {
            y
        }
    }

    /// Complementary error function
    ///
    /// Calculates the complementary probability.
    /// Is `1 - erf(x)`. Is computed directly, so that you can use it to avoid
    /// the loss of precision that would result from subtracting
    /// large probabilities (on large `x`) from 1.
    pub fn erfc(self) -> Self
    {
        /*
        * Coefficients for approximation to  erf on [0,0.84375]
        */
        const PP: [f64; 5] = [
            1.28379167095512558561e-01, /* 0x3FC06EBA, 0x8214DB68 */
            -3.25042107247001499370e-01, /* 0xBFD4CD7D, 0x691CB913 */
            -2.84817495755985104766e-02, /* 0xBF9D2A51, 0xDBD7194F */
            -5.77027029648944159157e-03, /* 0xBF77A291, 0x236668E4 */
            -2.37630166566501626084e-05 /* 0xBEF8EAD6, 0x120016AC */
        ];
        const QQ: [f64; 5] = [
            3.97917223959155352819e-01, /* 0x3FD97779, 0xCDDADC09 */
            6.50222499887672944485e-02, /* 0x3FB0A54C, 0x5536CEBA */
            5.08130628187576562776e-03, /* 0x3F74D022, 0xC4D36B0F */
            1.32494738004321644526e-04, /* 0x3F215DC9, 0x221C1A10 */
            -3.96022827877536812320e-06 /* 0xBED09C43, 0x42A26120 */
        ];

        let sign = self.is_sign_negative();
        if self.is_nan()
        {
            return self
        }
        let two = Self::from_uint(2u8);
        if self.is_infinite()
        {
            return if sign {two} else {Self::zero()}
        }
        let xabs = self.abs();
        let one = Self::one();
        if xabs < Self::from(0.84375)
        {
            /* |x| < 0.84375 */
            if xabs < Self::from(1.3877787807814456755295395851135e-17)
            {
                /* |x| < 2**-56 */
                return one - self;
            }
            let z = self*self;
            let mut b = false;
            let mut r = util::polynomial(&PP, z, SIGN_BIT, false);
            let mut s = one + z*util::polynomial(&QQ, z, SIGN_BIT, false);
            if !SIGN_BIT && r.is_nan()
            {
                r = util::polynomial(&PP.map(Neg::neg), z, SIGN_BIT, false);
                b = !b;
            }
            if !SIGN_BIT && s.is_nan()
            {
                s = z*util::polynomial(&QQ.map(Neg::neg), z, SIGN_BIT, false) - one;
                b = !b;
            }
            let y = r / s;
            let half = Self::from(0.5);
            if sign || xabs < Self::from(1.0/4.0) || (!SIGN_BIT && self < half)
            {
                /* x < 1/4 */
                if b
                {
                    if self < self*y
                    {
                        return one + (self*y - self);
                    }
                    return one - (self - self*y);
                }
                return one - (self + self*y);
            }
            if b
            {
                if self - half < self*y
                {
                    return half + (self*y + half - self);
                }
                return half - (self - half - self*y);
            }
            return half - (self - half + self*y);
        }
        if xabs < Self::from_uint(28u8)
        {
            /* 0.84375 <= |x| < 28 */
            if sign
            {
                return two - self.erfc2();
            }

            return self.erfc2();
        }

        let x1p_1022 = Self::from(f64::from_bits(0x0010000000000000));
        if sign
        {
            two - x1p_1022
        }
        else
        {
            x1p_1022 * x1p_1022
        }
    }
    
    fn bessel0_p(self) -> (Self, bool)
    {
        /* The asymptotic expansions of pzero is
        *      1 - 9/128 s^2 + 11025/98304 s^4 - ...,  where s = 1/x.
        * For x >= 2, We approximate pzero by
        *      pzero(x) = 1 + (R/S)
        * where  R = pR0 + pR1*s^2 + pR2*s^4 + ... + pR5*s^10
        *        S = 1 + pS0*s^2 + ... + pS4*s^10
        * and
        *      | pzero(x)-1-R/S | <= 2  ** ( -60.26)
        */
        const PR8: [f64; 6] = [
            /* for x in [inf, 8]=1/[0,0.125] */
            0.00000000000000000000e+00,  /* 0x00000000, 0x00000000 */
            -7.03124999999900357484e-02, /* 0xBFB1FFFF, 0xFFFFFD32 */
            -8.08167041275349795626e+00, /* 0xC02029D0, 0xB44FA779 */
            -2.57063105679704847262e+02, /* 0xC0701102, 0x7B19E863 */
            -2.48521641009428822144e+03, /* 0xC0A36A6E, 0xCD4DCAFC */
            -5.25304380490729545272e+03, /* 0xC0B4850B, 0x36CC643D */
        ];
        const PS8: [f64; 5] = [
            1.16534364619668181717e+02, /* 0x405D2233, 0x07A96751 */
            3.83374475364121826715e+03, /* 0x40ADF37D, 0x50596938 */
            4.05978572648472545552e+04, /* 0x40E3D2BB, 0x6EB6B05F */
            1.16752972564375915681e+05, /* 0x40FC810F, 0x8F9FA9BD */
            4.76277284146730962675e+04, /* 0x40E74177, 0x4F2C49DC */
        ];

        const PR5: [f64; 6] = [
            /* for x in [8,4.5454]=1/[0.125,0.22001] */
            -1.14125464691894502584e-11, /* 0xBDA918B1, 0x47E495CC */
            -7.03124940873599280078e-02, /* 0xBFB1FFFF, 0xE69AFBC6 */
            -4.15961064470587782438e+00, /* 0xC010A370, 0xF90C6BBF */
            -6.76747652265167261021e+01, /* 0xC050EB2F, 0x5A7D1783 */
            -3.31231299649172967747e+02, /* 0xC074B3B3, 0x6742CC63 */
            -3.46433388365604912451e+02, /* 0xC075A6EF, 0x28A38BD7 */
        ];
        const PS5: [f64; 5] = [
            6.07539382692300335975e+01, /* 0x404E6081, 0x0C98C5DE */
            1.05125230595704579173e+03, /* 0x40906D02, 0x5C7E2864 */
            5.97897094333855784498e+03, /* 0x40B75AF8, 0x8FBE1D60 */
            9.62544514357774460223e+03, /* 0x40C2CCB8, 0xFA76FA38 */
            2.40605815922939109441e+03, /* 0x40A2CC1D, 0xC70BE864 */
        ];

        const PR3: [f64; 6] = [
            /* for x in [4.547,2.8571]=1/[0.2199,0.35001] */
            -2.54704601771951915620e-09, /* 0xBE25E103, 0x6FE1AA86 */
            -7.03119616381481654654e-02, /* 0xBFB1FFF6, 0xF7C0E24B */
            -2.40903221549529611423e+00, /* 0xC00345B2, 0xAEA48074 */
            -2.19659774734883086467e+01, /* 0xC035F74A, 0x4CB94E14 */
            -5.80791704701737572236e+01, /* 0xC04D0A22, 0x420A1A45 */
            -3.14479470594888503854e+01, /* 0xC03F72AC, 0xA892D80F */
        ];
        const PS3: [f64; 5] = [
            3.58560338055209726349e+01, /* 0x4041ED92, 0x84077DD3 */
            3.61513983050303863820e+02, /* 0x40769839, 0x464A7C0E */
            1.19360783792111533330e+03, /* 0x4092A66E, 0x6D1061D6 */
            1.12799679856907414432e+03, /* 0x40919FFC, 0xB8C39B7E */
            1.73580930813335754692e+02, /* 0x4065B296, 0xFC379081 */
        ];

        const PR2: [f64; 6] = [
            /* for x in [2.8570,2]=1/[0.3499,0.5] */
            -8.87534333032526411254e-08, /* 0xBE77D316, 0xE927026D */
            -7.03030995483624743247e-02, /* 0xBFB1FF62, 0x495E1E42 */
            -1.45073846780952986357e+00, /* 0xBFF73639, 0x8A24A843 */
            -7.63569613823527770791e+00, /* 0xC01E8AF3, 0xEDAFA7F3 */
            -1.11931668860356747786e+01, /* 0xC02662E6, 0xC5246303 */
            -3.23364579351335335033e+00, /* 0xC009DE81, 0xAF8FE70F */
        ];
        const PS2: [f64; 5] = [
            2.22202997532088808441e+01, /* 0x40363865, 0x908B5959 */
            1.36206794218215208048e+02, /* 0x4061069E, 0x0EE8878F */
            2.70470278658083486789e+02, /* 0x4070E786, 0x42EA079B */
            1.53875394208320329881e+02, /* 0x40633C03, 0x3AB6FAFF */
            1.46576176948256193810e+01, /* 0x402D50B3, 0x44391809 */
        ];

        let p: &[f64; 6];
        let q: &[f64; 5];

        let xabs = self.abs();
        if xabs >= Self::from_uint(8u8)
        {
            p = &PR8;
            q = &PS8;
        }
        else if xabs >= Self::from(4.5454)
        {
            p = &PR5;
            q = &PS5;
        }
        else if xabs >= Self::from(2.857)
        {
            p = &PR3;
            q = &PS3;
        }
        else
        {
            p = &PR2;
            q = &PS2;
        }
        let one = Self::one();
        let z = (self*self).recip();
        let mut b = false;
        let mut r = util::polynomial(p, z, SIGN_BIT, false);
        let mut s = one + z*util::polynomial(q, z, SIGN_BIT, false);
        if !SIGN_BIT && r.is_nan()
        {
            b = !b;
            r = util::polynomial(&p.map(Neg::neg), z, SIGN_BIT, false);
        }
        if !SIGN_BIT && s.is_nan()
        {
            s = z*util::polynomial(&q.map(Neg::neg), z, SIGN_BIT, false);
            if s >= one
            {
                b = !b;
                s -= one
            }
            else
            {
                s = one - s
            }
        }
        
        let rs = r/s;
        if b
        {
            if one < rs
            {
                return (rs - one, true)
            }
            return (one - rs, false)
        }
        (one + rs, false)
    }
    
    fn bessel0_q(self) -> (Self, bool)
    {
        /* For x >= 8, the asymptotic expansions of qzero is
        *      -1/8 s + 75/1024 s^3 - ..., where s = 1/x.
        * We approximate pzero by
        *      qzero(x) = s*(-1.25 + (R/S))
        * where  R = qR0 + qR1*s^2 + qR2*s^4 + ... + qR5*s^10
        *        S = 1 + qS0*s^2 + ... + qS5*s^12
        * and
        *      | qzero(x)/s +1.25-R/S | <= 2  ** ( -61.22)
        */
        const QR8: [f64; 6] = [
            /* for x in [inf, 8]=1/[0,0.125] */
            0.00000000000000000000e+00, /* 0x00000000, 0x00000000 */
            7.32421874999935051953e-02, /* 0x3FB2BFFF, 0xFFFFFE2C */
            1.17682064682252693899e+01, /* 0x40278952, 0x5BB334D6 */
            5.57673380256401856059e+02, /* 0x40816D63, 0x15301825 */
            8.85919720756468632317e+03, /* 0x40C14D99, 0x3E18F46D */
            3.70146267776887834771e+04, /* 0x40E212D4, 0x0E901566 */
        ];
        const QS8: [f64; 6] = [
            1.63776026895689824414e+02,  /* 0x406478D5, 0x365B39BC */
            8.09834494656449805916e+03,  /* 0x40BFA258, 0x4E6B0563 */
            1.42538291419120476348e+05,  /* 0x41016652, 0x54D38C3F */
            8.03309257119514397345e+05,  /* 0x412883DA, 0x83A52B43 */
            8.40501579819060512818e+05,  /* 0x4129A66B, 0x28DE0B3D */
            -3.43899293537866615225e+05, /* 0xC114FD6D, 0x2C9530C5 */
        ];

        const QR5: [f64; 6] = [
            /* for x in [8,4.5454]=1/[0.125,0.22001] */
            1.84085963594515531381e-11, /* 0x3DB43D8F, 0x29CC8CD9 */
            7.32421766612684765896e-02, /* 0x3FB2BFFF, 0xD172B04C */
            5.83563508962056953777e+00, /* 0x401757B0, 0xB9953DD3 */
            1.35111577286449829671e+02, /* 0x4060E392, 0x0A8788E9 */
            1.02724376596164097464e+03, /* 0x40900CF9, 0x9DC8C481 */
            1.98997785864605384631e+03, /* 0x409F17E9, 0x53C6E3A6 */
        ];
        const QS5: [f64; 6] = [
            8.27766102236537761883e+01,  /* 0x4054B1B3, 0xFB5E1543 */
            2.07781416421392987104e+03,  /* 0x40A03BA0, 0xDA21C0CE */
            1.88472887785718085070e+04,  /* 0x40D267D2, 0x7B591E6D */
            5.67511122894947329769e+04,  /* 0x40EBB5E3, 0x97E02372 */
            3.59767538425114471465e+04,  /* 0x40E19118, 0x1F7A54A0 */
            -5.35434275601944773371e+03, /* 0xC0B4EA57, 0xBEDBC609 */
        ];

        const QR3: [f64; 6] = [
            /* for x in [4.547,2.8571]=1/[0.2199,0.35001] */
            4.37741014089738620906e-09, /* 0x3E32CD03, 0x6ADECB82 */
            7.32411180042911447163e-02, /* 0x3FB2BFEE, 0x0E8D0842 */
            3.34423137516170720929e+00, /* 0x400AC0FC, 0x61149CF5 */
            4.26218440745412650017e+01, /* 0x40454F98, 0x962DAEDD */
            1.70808091340565596283e+02, /* 0x406559DB, 0xE25EFD1F */
            1.66733948696651168575e+02, /* 0x4064D77C, 0x81FA21E0 */
        ];
        const QS3: [f64; 6] = [
            4.87588729724587182091e+01,  /* 0x40486122, 0xBFE343A6 */
            7.09689221056606015736e+02,  /* 0x40862D83, 0x86544EB3 */
            3.70414822620111362994e+03,  /* 0x40ACF04B, 0xE44DFC63 */
            6.46042516752568917582e+03,  /* 0x40B93C6C, 0xD7C76A28 */
            2.51633368920368957333e+03,  /* 0x40A3A8AA, 0xD94FB1C0 */
            -1.49247451836156386662e+02, /* 0xC062A7EB, 0x201CF40F */
        ];

        const QR2: [f64; 6] = [
            /* for x in [2.8570,2]=1/[0.3499,0.5] */
            1.50444444886983272379e-07, /* 0x3E84313B, 0x54F76BDB */
            7.32234265963079278272e-02, /* 0x3FB2BEC5, 0x3E883E34 */
            1.99819174093815998816e+00, /* 0x3FFFF897, 0xE727779C */
            1.44956029347885735348e+01, /* 0x402CFDBF, 0xAAF96FE5 */
            3.16662317504781540833e+01, /* 0x403FAA8E, 0x29FBDC4A */
            1.62527075710929267416e+01, /* 0x403040B1, 0x71814BB4 */
        ];
        const QS2: [f64; 6] = [
            3.03655848355219184498e+01,  /* 0x403E5D96, 0xF7C07AED */
            2.69348118608049844624e+02,  /* 0x4070D591, 0xE4D14B40 */
            8.44783757595320139444e+02,  /* 0x408A6645, 0x22B3BF22 */
            8.82935845112488550512e+02,  /* 0x408B977C, 0x9C5CC214 */
            2.12666388511798828631e+02,  /* 0x406A9553, 0x0E001365 */
            -5.31095493882666946917e+00, /* 0xC0153E6A, 0xF8B32931 */
        ];

        let p: &[f64; 6];
        let q: &[f64; 6];

        let xabs = self.abs();
        if xabs >= Self::from_uint(8u8)
        {
            p = &QR8;
            q = &QS8;
        }
        else if xabs >= Self::from(4.5454)
        {
            p = &QR5;
            q = &QS5;
        }
        else if xabs >= Self::from(2.857)
        {
            p = &QR3;
            q = &QS3;
        }
        else
        {
            p = &QR2;
            q = &QS2;
        }
        let one = Self::one();
        let z = (self*self).recip();
        let mut b = false;
        let mut r = util::polynomial(p, z, SIGN_BIT, false);
        let mut s = one + z*util::polynomial(q, z, SIGN_BIT, false);
        if !SIGN_BIT && r.is_nan()
        {
            b = !b;
            r = util::polynomial(&p.map(Neg::neg), z, SIGN_BIT, false);
        }
        if !SIGN_BIT && s.is_nan()
        {
            b = !b;
            s = z*util::polynomial(&q.map(Neg::neg), z, SIGN_BIT, false) - one;
        }

        let rs = r/s;
        let ofs = Self::from(0.125);
        if b
        {
            return ((rs + ofs)/self, true)
        }
        
        if !SIGN_BIT && rs < ofs
        {
            return ((ofs - rs)/self, true)
        }
        ((rs - ofs)/self, false)
    }

    fn bessel0_common(self, y0: bool) -> (Self, bool)
    {
        const INVSQRTPI: f64 = 5.64189583547756279280e-01;

        /*
         * j0(x) = sqrt(2/(pi*x))*(p0(x)*cos(x-pi/4)-q0(x)*sin(x-pi/4))
         * y0(x) = sqrt(2/(pi*x))*(p0(x)*sin(x-pi/4)+q0(x)*cos(x-pi/4))
         *
         * sin(x-pi/4) = (sin(x) - cos(x))/sqrt(2)
         * cos(x-pi/4) = (sin(x) + cos(x))/sqrt(2)
         * sin(x) +- cos(x) = -cos(2x)/(sin(x) -+ cos(x))
         */
        let (s, s_s) = self.sin_extra_sign();
        let (mut c, mut c_s) = self.cos_extra_sign();
        if y0
        {
            if SIGN_BIT
            {
                c = -c;
            }
            else
            {
                c_s = !c_s;
            }
        }
        
        let mut cc_s = false;
        let mut cc = match(s_s, c_s)
        {
            (false, false) => s + c,
            (true, false) => if SIGN_BIT || s <= c
            {
                c - s
            }
            else
            {
                cc_s = !cc_s;
                s - c
            },
            (false, true) => if SIGN_BIT || c <= s
            {
                s - c
            }
            else
            {
                cc_s = !cc_s;
                c - s
            },
            (true, true) => if SIGN_BIT
            {
                -s - c
            }
            else
            {
                cc_s = !cc_s;
                s + c
            }
        };

        let two = Self::from_uint(2u8);
        let xabs = self.abs();
        /* avoid overflow in 2*x, big ulp error when x>=0x1p1023 */
        if xabs < Self::max_value()/two
        {
            let mut ss_s = false;
            let mut ss = match(s_s, !c_s)
            {
                (false, false) => s + c,
                (true, false) => if SIGN_BIT || s <= c
                {
                    c - s
                }
                else
                {
                    ss_s = !ss_s;
                    s - c
                },
                (false, true) => if SIGN_BIT || c <= s
                {
                    s - c
                }
                else
                {
                    ss_s = !ss_s;
                    c - s
                },
                (true, true) => if SIGN_BIT
                {
                    -s - c
                }
                else
                {
                    ss_s = !ss_s;
                    s + c
                }
            };

            let (mut z, mut z_s) = (two*self).cos_extra_sign();
            if SIGN_BIT
            {
                z = -z
            }
            else
            {
                z_s = !z_s
            }

            if if s_s == c_s
            {
                s*c < Self::zero()
            }
            else
            {
                s*c > Self::zero()
            }
            {
                cc_s = z_s^ss_s;
                cc = z / ss;
            }
            else
            {
                ss_s = z_s^cc_s;
                ss = z / cc;
            }
            if xabs < Self::from((EXP_BASE as f64).powf(17.0/127.0*Self::exp_bias().to_f64().unwrap()))
            {
                if y0
                {
                    if SIGN_BIT
                    {
                        ss = -ss;
                    }
                    else
                    {
                        ss_s = !ss_s;
                    }
                }
                
                let (mut p, mut p_s) = self.bessel0_p();
                p *= cc;
                p_s ^= cc_s;

                let (mut q, mut q_s) = self.bessel0_q();
                q *= ss;
                q_s ^= ss_s;

                cc_s = false;
                cc = match(p_s, !q_s)
                {
                    (false, false) => p + q,
                    (true, false) => if SIGN_BIT || p <= q
                    {
                        q - p
                    }
                    else
                    {
                        cc_s = !cc_s;
                        p - q
                    },
                    (false, true) => if SIGN_BIT || q <= p
                    {
                        p - q
                    }
                    else
                    {
                        cc_s = !cc_s;
                        q - p
                    },
                    (true, true) => if SIGN_BIT
                    {
                        -p - q
                    }
                    else
                    {
                        cc_s = !cc_s;
                        p + q
                    }
                };
            }
        }
        
        (Self::from(INVSQRTPI)*cc / self.sqrt(), cc_s)
    }

    fn j0_extra_sign(mut self) -> (Self, bool)
    {
        self = self.abs();
        /* j0(+-inf)=0, j0(nan)=nan */
        if !self.is_finite()
        {
            return ((self*self).recip(), false)
        }

        let two = Self::from_uint(2u8);
        if self >= two
        {
            /* |x| >= 2 */
            /* large ulp error near zeros: 2.4, 5.52, 8.6537,.. */
            return self.bessel0_common(false)
        }

        let one = Self::one();

        /* 1 - x*x/4 + x*x*R(x^2)/S(x^2) */
        if self >= Self::from(0.0001220703125)
        {
            /* R0/S0 on [0, 2.00] */
            const R0: [f64; 4] = [
                1.56249999999999947958e-02, /* 0x3F8FFFFF, 0xFFFFFFFD */
                -1.89979294238854721751e-04, /* 0xBF28E6A5, 0xB61AC6E9 */
                1.82954049532700665670e-06, /* 0x3EBEB1D1, 0x0C503919 */
                -4.61832688532103189199e-09 /* 0xBE33D5E7, 0x73D63FCE */
            ];
            const S0: [f64; 4] = [
                1.56191029464890010492e-02, /* 0x3F8FFCE8, 0x82C8C2A4 */
                1.16926784663337450260e-04, /* 0x3F1EA6D2, 0xDD57DBF4 */
                5.13546550207318111446e-07, /* 0x3EA13B54, 0xCE84D5A9 */
                1.16614003333790000205e-09 /* 0x3E1408BC, 0xF4745D8F */
            ];

            /* |x| >= 2**-13 */
            /* up to 4ulp error close to 2 */
            let z = self*self;
            let mut b = false;
            let mut r = z*util::polynomial(&R0, z, SIGN_BIT, false);
            let mut s = one + z*util::polynomial(&S0, z, SIGN_BIT, false);
            if !SIGN_BIT && r.is_nan()
            {
                b = !b;
                r = z*util::polynomial(&R0.map(Neg::neg), z, SIGN_BIT, false)
            }
            if !SIGN_BIT && s.is_nan()
            {
                s = z*util::polynomial(&S0.map(Neg::neg), z, SIGN_BIT, false);
                if s > one
                {
                    b = !b;
                    s -= one
                }
                else
                {
                    s = one - s
                }
            }

            let rs = r/s;
            let zrs = z*rs;
            let d = (one + self/two)*(one - self/two);

            if b
            {
                if d < zrs
                {
                    return (zrs - d, true)
                }
                return (d - zrs, false);
            }
            return (d + zrs, false);
        }

        /* 1 - x*x/4 */
        /* prevent underflow */
        /* inexact should be raised when x!=0, this is not done correctly */
        if self >= Self::from((EXP_BASE as f64).powf(-13.0/127.0*Self::exp_bias().to_f64().unwrap()))
        {
            /* |x| >= 2**-127 */
            self = Self::from(0.25)*self*self;
        }

        if one < self
        {
            return (self - one, true)
        }
        (one - self, false)
    }
    
    /// Bessel function of the first kind with α = 0.
    pub fn j0(self) -> Self
    {
        let (mut y, s) = self.j0_extra_sign();
        if s
        {
            y = -y
        }
        y
    }

    fn y0_extra_sign(self) -> (Self, bool)
    {
        /* y0(nan)=nan, y0(<0)=nan, y0(0)=-inf, y0(inf)=0 */
        if self.is_zero()
        {
            if !SIGN_BIT
            {
                return (Self::infinity(), true)
            }
            return (Self::neg_infinity(), false)
        }
        if self.is_sign_negative()
        {
            return (Self::snan(), false)
        }
        if !self.is_finite()
        {
            return (self.recip(), false)
        }

        if self >= Self::from_uint(2u8)
        {
            /* x >= 2 */
            /* large ulp errors near zeros: 3.958, 7.086,.. */
            return self.bessel0_common(true);
        }

        const TPI: f64 = 6.36619772367581382433e-01; /* 0x3FE45F30, 0x6DC9C883 */
        
        const U0: [f64; 7] = [
            -7.38042951086872317523e-02, /* 0xBFB2E4D6, 0x99CBD01F */
            1.76666452509181115538e-01, /* 0x3FC69D01, 0x9DE9E3FC */
            -1.38185671945596898896e-02, /* 0xBF8C4CE8, 0xB16CFA97 */
            3.47453432093683650238e-04, /* 0x3F36C54D, 0x20B29B6B */
            -3.81407053724364161125e-06, /* 0xBECFFEA7, 0x73D25CAD */
            1.95590137035022920206e-08, /* 0x3E550057, 0x3B4EABD4 */
            -3.98205194132103398453e-11 /* 0xBDC5E43D, 0x693FB3C8 */
        ];
        const V0: [f64; 4] = [
            1.27304834834123699328e-02, /* 0x3F8A1270, 0x91C9C71A */
            7.60068627350353253702e-05, /* 0x3F13ECBB, 0xF578C6C1 */
            2.59150851840457805467e-07, /* 0x3E91642D, 0x7FF202FD */
            4.41110311332675467403e-10 /* 0x3DFE5018, 0x3BD6D9EF */
        ];

        /* U(x^2)/V(x^2) + (2/pi)*j0(x)*log(x) */
        if self >= Self::from((EXP_BASE as f64).powf(-1.625*EXP_SIZE as f64))
        {
            /* large ulp error near the first zero, x ~= 0.89 */
            let z = self*self;

            let mut b = false;
            let mut u = util::polynomial(&U0, z, SIGN_BIT, false);
            let mut v = Self::one() + z*util::polynomial(&V0, z, SIGN_BIT, false);
            if !SIGN_BIT && u.is_nan()
            {
                b = !b;
                u = util::polynomial(&U0.map(Neg::neg), z, SIGN_BIT, false);
            }
            let one = Self::one();
            if !SIGN_BIT && v.is_nan()
            {
                v = z*util::polynomial(&V0.map(Neg::neg), z, SIGN_BIT, false);
                if v > one
                {
                    b = !b;
                    v -= one;
                }
                else
                {
                    v = one - v;
                }
            }

            let uv = u/v;
            if b
            {
                let (xj0, xj0_s) = self.j0_extra_sign();
                if !SIGN_BIT && self < one
                {
                    let l = Self::from(TPI)*(xj0*self.recip().ln());
                    if xj0_s
                    {
                        if l < uv
                        {
                            return (uv - l, true)
                        }
                        return (l - uv, false);
                    }
                    return (uv + l, true);
                }
                let l = Self::from(TPI)*(xj0*self.ln());

                if xj0_s
                {
                    return (l + uv, true);
                }
                if l < uv
                {
                    return (uv - l, true)
                }
                return (l - uv, false);
            }
            if !SIGN_BIT && self < one
            {
                let (xj0, xj0_s) = self.j0_extra_sign();
                let l = Self::from(TPI)*(xj0*self.recip().ln());
                
                if xj0_s
                {
                    return (uv + l, false);
                }
                if uv < l
                {
                    return (l - uv, true)
                }
                return (uv - l, false);
            }
            let (xj0, xj0_s) = self.j0_extra_sign();
            let l = Self::from(TPI)*(xj0*self.ln());
            if xj0_s
            {
                if uv < l
                {
                    return (l - uv, true)
                }
                return (uv - l, false);
            }
            return (uv + l, false);
        }

        let u00 = Self::from(U0[0]);
        if !SIGN_BIT
        {
            let l = Self::from(TPI)*self.recip().ln();

            if u00 < l
            {
                return (l - u00, true)
            }
            return (u00 - l, false)
        }
        (u00 + Self::from(TPI)*self.ln(), false)
    }

    /// Bessel function of the second kind with α = 0.
    pub fn y0(self) -> Self
    {
        let (mut y, s) = self.y0_extra_sign();
        if s
        {
            y = -y
        }
        y
    }
    
    fn bessel1_p(self) -> (Self, bool)
    {
        /* For x >= 8, the asymptotic expansions of pone is
        *      1 + 15/128 s^2 - 4725/2^15 s^4 - ...,   where s = 1/x.
        * We approximate pone by
        *      pone(x) = 1 + (R/S)
        * where  R = pr0 + pr1*s^2 + pr2*s^4 + ... + pr5*s^10
        *        S = 1 + ps0*s^2 + ... + ps4*s^10
        * and
        *      | pone(x)-1-R/S | <= 2  ** ( -60.06)
        */

        const PR8: [f64; 6] = [
            /* for x in [inf, 8]=1/[0,0.125] */
            0.00000000000000000000e+00, /* 0x00000000, 0x00000000 */
            1.17187499999988647970e-01, /* 0x3FBDFFFF, 0xFFFFFCCE */
            1.32394806593073575129e+01, /* 0x402A7A9D, 0x357F7FCE */
            4.12051854307378562225e+02, /* 0x4079C0D4, 0x652EA590 */
            3.87474538913960532227e+03, /* 0x40AE457D, 0xA3A532CC */
            7.91447954031891731574e+03, /* 0x40BEEA7A, 0xC32782DD */
        ];
        const PS8: [f64; 5] = [
            1.14207370375678408436e+02, /* 0x405C8D45, 0x8E656CAC */
            3.65093083420853463394e+03, /* 0x40AC85DC, 0x964D274F */
            3.69562060269033463555e+04, /* 0x40E20B86, 0x97C5BB7F */
            9.76027935934950801311e+04, /* 0x40F7D42C, 0xB28F17BB */
            3.08042720627888811578e+04, /* 0x40DE1511, 0x697A0B2D */
        ];

        const PR5: [f64; 6] = [
            /* for x in [8,4.5454]=1/[0.125,0.22001] */
            1.31990519556243522749e-11, /* 0x3DAD0667, 0xDAE1CA7D */
            1.17187493190614097638e-01, /* 0x3FBDFFFF, 0xE2C10043 */
            6.80275127868432871736e+00, /* 0x401B3604, 0x6E6315E3 */
            1.08308182990189109773e+02, /* 0x405B13B9, 0x452602ED */
            5.17636139533199752805e+02, /* 0x40802D16, 0xD052D649 */
            5.28715201363337541807e+02, /* 0x408085B8, 0xBB7E0CB7 */
        ];
        const PS5: [f64; 5] = [
            5.92805987221131331921e+01, /* 0x404DA3EA, 0xA8AF633D */
            9.91401418733614377743e+02, /* 0x408EFB36, 0x1B066701 */
            5.35326695291487976647e+03, /* 0x40B4E944, 0x5706B6FB */
            7.84469031749551231769e+03, /* 0x40BEA4B0, 0xB8A5BB15 */
            1.50404688810361062679e+03, /* 0x40978030, 0x036F5E51 */
        ];

        const PR3: [f64; 6] = [
            3.02503916137373618024e-09, /* 0x3E29FC21, 0xA7AD9EDD */
            1.17186865567253592491e-01, /* 0x3FBDFFF5, 0x5B21D17B */
            3.93297750033315640650e+00, /* 0x400F76BC, 0xE85EAD8A */
            3.51194035591636932736e+01, /* 0x40418F48, 0x9DA6D129 */
            9.10550110750781271918e+01, /* 0x4056C385, 0x4D2C1837 */
            4.85590685197364919645e+01, /* 0x4048478F, 0x8EA83EE5 */
        ];
        const PS3: [f64; 5] = [
            3.47913095001251519989e+01, /* 0x40416549, 0xA134069C */
            3.36762458747825746741e+02, /* 0x40750C33, 0x07F1A75F */
            1.04687139975775130551e+03, /* 0x40905B7C, 0x5037D523 */
            8.90811346398256432622e+02, /* 0x408BD67D, 0xA32E31E9 */
            1.03787932439639277504e+02, /* 0x4059F26D, 0x7C2EED53 */
        ];

        const PR2: [f64; 6] = [
            /* for x in [2.8570,2]=1/[0.3499,0.5] */
            1.07710830106873743082e-07, /* 0x3E7CE9D4, 0xF65544F4 */
            1.17176219462683348094e-01, /* 0x3FBDFF42, 0xBE760D83 */
            2.36851496667608785174e+00, /* 0x4002F2B7, 0xF98FAEC0 */
            1.22426109148261232917e+01, /* 0x40287C37, 0x7F71A964 */
            1.76939711271687727390e+01, /* 0x4031B1A8, 0x177F8EE2 */
            5.07352312588818499250e+00, /* 0x40144B49, 0xA574C1FE */
        ];
        const PS2: [f64; 5] = [
            2.14364859363821409488e+01, /* 0x40356FBD, 0x8AD5ECDC */
            1.25290227168402751090e+02, /* 0x405F5293, 0x14F92CD5 */
            2.32276469057162813669e+02, /* 0x406D08D8, 0xD5A2DBD9 */
            1.17679373287147100768e+02, /* 0x405D6B7A, 0xDA1884A9 */
            8.36463893371618283368e+00, /* 0x4020BAB1, 0xF44E5192 */
        ];

        let p: &[f64; 6];
        let q: &[f64; 5];
        
        let xabs = self.abs();
        if xabs >= Self::from_uint(8u8)
        {
            p = &PR8;
            q = &PS8;
        }
        else if xabs >= Self::from(4.5454)
        {
            p = &PR5;
            q = &PS5;
        }
        else if xabs >= Self::from(2.857)
        {
            p = &PR3;
            q = &PS3;
        }
        else
        {
            p = &PR2;
            q = &PS2;
        }
        let one = Self::one();
        let z = (self*self).recip();

        let mut b = false;
        let mut r = util::polynomial(p, z, SIGN_BIT, false);
        let mut s = one + z*util::polynomial(q, z, SIGN_BIT, false);
        if !SIGN_BIT && r.is_nan()
        {
            b = !b;
            r = util::polynomial(&p.map(Neg::neg), z, SIGN_BIT, false);
        }
        if !SIGN_BIT && s.is_nan()
        {
            s = z*util::polynomial(&q.map(Neg::neg), z, SIGN_BIT, false);
            if s > one
            {
                b = !b;
                s -= one
            }
            else
            {
                s = one - s
            }
        }
        
        let rs = r/s;
        if b
        {
            if rs > one
            {
                return (rs - one, true)
            }
            return (one - rs, false)
        }
        (one + rs, false)
    }

    fn bessel1_q(self) -> (Self, bool)
    {
        /* For x >= 8, the asymptotic expansions of qone is
        *      3/8 s - 105/1024 s^3 - ..., where s = 1/x.
        * We approximate pone by
        *      qone(x) = s*(0.375 + (R/S))
        * where  R = qr1*s^2 + qr2*s^4 + ... + qr5*s^10
        *        S = 1 + qs1*s^2 + ... + qs6*s^12
        * and
        *      | qone(x)/s -0.375-R/S | <= 2  ** ( -61.13)
        */

        const QR8: [f64; 6] = [
            /* for x in [inf, 8]=1/[0,0.125] */
            0.00000000000000000000e+00,  /* 0x00000000, 0x00000000 */
            -1.02539062499992714161e-01, /* 0xBFBA3FFF, 0xFFFFFDF3 */
            -1.62717534544589987888e+01, /* 0xC0304591, 0xA26779F7 */
            -7.59601722513950107896e+02, /* 0xC087BCD0, 0x53E4B576 */
            -1.18498066702429587167e+04, /* 0xC0C724E7, 0x40F87415 */
            -4.84385124285750353010e+04, /* 0xC0E7A6D0, 0x65D09C6A */
        ];
        const QS8: [f64; 6] = [
            1.61395369700722909556e+02,  /* 0x40642CA6, 0xDE5BCDE5 */
            7.82538599923348465381e+03,  /* 0x40BE9162, 0xD0D88419 */
            1.33875336287249578163e+05,  /* 0x4100579A, 0xB0B75E98 */
            7.19657723683240939863e+05,  /* 0x4125F653, 0x72869C19 */
            6.66601232617776375264e+05,  /* 0x412457D2, 0x7719AD5C */
            -2.94490264303834643215e+05, /* 0xC111F969, 0x0EA5AA18 */
        ];

        const QR5: [f64; 6] = [
            /* for x in [8,4.5454]=1/[0.125,0.22001] */
            -2.08979931141764104297e-11, /* 0xBDB6FA43, 0x1AA1A098 */
            -1.02539050241375426231e-01, /* 0xBFBA3FFF, 0xCB597FEF */
            -8.05644828123936029840e+00, /* 0xC0201CE6, 0xCA03AD4B */
            -1.83669607474888380239e+02, /* 0xC066F56D, 0x6CA7B9B0 */
            -1.37319376065508163265e+03, /* 0xC09574C6, 0x6931734F */
            -2.61244440453215656817e+03, /* 0xC0A468E3, 0x88FDA79D */
        ];
        const QS5: [f64; 6] = [
            8.12765501384335777857e+01,  /* 0x405451B2, 0xFF5A11B2 */
            1.99179873460485964642e+03,  /* 0x409F1F31, 0xE77BF839 */
            1.74684851924908907677e+04,  /* 0x40D10F1F, 0x0D64CE29 */
            4.98514270910352279316e+04,  /* 0x40E8576D, 0xAABAD197 */
            2.79480751638918118260e+04,  /* 0x40DB4B04, 0xCF7C364B */
            -4.71918354795128470869e+03, /* 0xC0B26F2E, 0xFCFFA004 */
        ];

        const QR3: [f64; 6] = [
            -5.07831226461766561369e-09, /* 0xBE35CFA9, 0xD38FC84F */
            -1.02537829820837089745e-01, /* 0xBFBA3FEB, 0x51AEED54 */
            -4.61011581139473403113e+00, /* 0xC01270C2, 0x3302D9FF */
            -5.78472216562783643212e+01, /* 0xC04CEC71, 0xC25D16DA */
            -2.28244540737631695038e+02, /* 0xC06C87D3, 0x4718D55F */
            -2.19210128478909325622e+02, /* 0xC06B66B9, 0x5F5C1BF6 */
        ];
        const QS3: [f64; 6] = [
            4.76651550323729509273e+01,  /* 0x4047D523, 0xCCD367E4 */
            6.73865112676699709482e+02,  /* 0x40850EEB, 0xC031EE3E */
            3.38015286679526343505e+03,  /* 0x40AA684E, 0x448E7C9A */
            5.54772909720722782367e+03,  /* 0x40B5ABBA, 0xA61D54A6 */
            1.90311919338810798763e+03,  /* 0x409DBC7A, 0x0DD4DF4B */
            -1.35201191444307340817e+02, /* 0xC060E670, 0x290A311F */
        ];

        const QR2: [f64; 6] = [
            /* for x in [2.8570,2]=1/[0.3499,0.5] */
            -1.78381727510958865572e-07, /* 0xBE87F126, 0x44C626D2 */
            -1.02517042607985553460e-01, /* 0xBFBA3E8E, 0x9148B010 */
            -2.75220568278187460720e+00, /* 0xC0060484, 0x69BB4EDA */
            -1.96636162643703720221e+01, /* 0xC033A9E2, 0xC168907F */
            -4.23253133372830490089e+01, /* 0xC04529A3, 0xDE104AAA */
            -2.13719211703704061733e+01, /* 0xC0355F36, 0x39CF6E52 */
        ];
        const QS2: [f64; 6] = [
            2.95333629060523854548e+01,  /* 0x403D888A, 0x78AE64FF */
            2.52981549982190529136e+02,  /* 0x406F9F68, 0xDB821CBA */
            7.57502834868645436472e+02,  /* 0x4087AC05, 0xCE49A0F7 */
            7.39393205320467245656e+02,  /* 0x40871B25, 0x48D4C029 */
            1.55949003336666123687e+02,  /* 0x40637E5E, 0x3C3ED8D4 */
            -4.95949898822628210127e+00, /* 0xC013D686, 0xE71BE86B */
        ];

        let p: &[f64; 6];
        let q: &[f64; 6];
        
        let xabs = self.abs();
        if xabs >= Self::from_uint(8u8)
        {
            p = &QR8;
            q = &QS8;
        }
        else if xabs >= Self::from(4.5454)
        {
            p = &QR5;
            q = &QS5;
        }
        else if xabs >= Self::from(2.857)
        {
            p = &QR3;
            q = &QS3;
        }
        else
        {
            p = &QR2;
            q = &QS2;
        }
        let one = Self::one();
        let z = (self*self).recip();
        let mut b = false;
        let mut r = util::polynomial(p, z, SIGN_BIT, false);
        let mut s = one + z*util::polynomial(q, z, SIGN_BIT, false);
        if !SIGN_BIT && r.is_nan()
        {
            b = !b;
            r = util::polynomial(&p.map(Neg::neg), z, SIGN_BIT, false);
        }
        if !SIGN_BIT && s.is_nan()
        {
            b = !b;
            s = z*util::polynomial(&q.map(Neg::neg), z, SIGN_BIT, false) - one;
        }

        let rs = r/s;
        let ofs = Self::from(0.375);
        if b
        {
            if rs > ofs
            {
                return ((rs - ofs)/self, true)
            }
            return ((ofs - rs)/self, false)
        }
        
        ((ofs + rs)/self, false)
    }
    
    fn bessel1_common(self, y1: bool, sign: bool) -> (Self, bool)
    {
        const INVSQRTPI: f64 = 5.64189583547756279280e-01;

        /*
        * j1(x) = sqrt(2/(pi*x))*(p1(x)*cos(x-3pi/4)-q1(x)*sin(x-3pi/4))
        * y1(x) = sqrt(2/(pi*x))*(p1(x)*sin(x-3pi/4)+q1(x)*cos(x-3pi/4))
        *
        * sin(x-3pi/4) = -(sin(x) + cos(x))/sqrt(2)
        * cos(x-3pi/4) = (sin(x) - cos(x))/sqrt(2)
        * sin(x) +- cos(x) = -cos(2x)/(sin(x) -+ cos(x))
        */
        let (mut s, mut s_s) = self.sin_extra_sign();
        if y1
        {
            if SIGN_BIT
            {
                s = -s;
            }
            else
            {
                s_s = !s_s;
            }
        }
        let (c, c_s) = self.cos_extra_sign();
        
        let mut cc_s = false;
        let mut cc = match(s_s, !c_s)
        {
            (false, false) => s + c,
            (true, false) => if SIGN_BIT || s <= c
            {
                c - s
            }
            else
            {
                cc_s = !cc_s;
                s - c
            },
            (false, true) => if SIGN_BIT || c <= s
            {
                s - c
            }
            else
            {
                cc_s = !cc_s;
                c - s
            },
            (true, true) => if SIGN_BIT
            {
                -s - c
            }
            else
            {
                cc_s = !cc_s;
                s + c
            }
        };

        let two = Self::from(2u8);
        let xabs = self.abs();
        if xabs < Self::max_value()/two
        {
            /* avoid overflow in 2*x */
            let mut ss_s = false;
            let mut ss = match(!s_s, !c_s)
            {
                (false, false) => s + c,
                (true, false) => if SIGN_BIT || s <= c
                {
                    c - s
                }
                else
                {
                    ss_s = !ss_s;
                    s - c
                },
                (false, true) => if SIGN_BIT || c <= s
                {
                    s - c
                }
                else
                {
                    ss_s = !ss_s;
                    c - s
                },
                (true, true) => if SIGN_BIT
                {
                    -s - c
                }
                else
                {
                    ss_s = !ss_s;
                    s + c
                }
            };
            
            let (z, z_s) = (two*self).cos_extra_sign();
            
            if if s_s == c_s
            {
                s*c < Self::zero()
            }
            else
            {
                s*c > Self::zero()
            }
            {
                cc_s = z_s^ss_s;
                cc = z / ss;
            }
            else
            {
                ss_s = z_s^cc_s;
                ss = z / cc;
            }

            if xabs < Self::from((EXP_BASE as f64).powf(17.0/127.0*Self::exp_bias().to_f64().unwrap()))
            {
                if y1
                {
                    if SIGN_BIT
                    {
                        ss = -ss;
                    }
                    else
                    {
                        ss_s = !ss_s;
                    }
                }
                
                let (mut p, mut p_s) = self.bessel1_p();
                p *= cc;
                p_s ^= cc_s;

                let (mut q, mut q_s) = self.bessel1_q();
                q *= ss;
                q_s ^= ss_s;

                cc_s = false;
                cc = match(p_s, !q_s)
                {
                    (false, false) => p + q,
                    (true, false) => if SIGN_BIT || p <= q
                    {
                        q - p
                    }
                    else
                    {
                        cc_s = !cc_s;
                        p - q
                    },
                    (false, true) => if SIGN_BIT || q <= p
                    {
                        p - q
                    }
                    else
                    {
                        cc_s = !cc_s;
                        q - p
                    },
                    (true, true) => if SIGN_BIT
                    {
                        -p - q
                    }
                    else
                    {
                        cc_s = !cc_s;
                        p + q
                    }
                };
            }
        }
        if sign
        {
            cc = -cc;
        }

        (Self::from(INVSQRTPI)*cc/self.sqrt(), cc_s)
    }

    fn j1_extra_sign(self) -> (Self, bool)
    {
        let sign = !self.sign_bit().is_zero();
        let xabs = self.abs();

        if xabs.is_infinite()
        {
            return ((self*self).recip(), false);
        }
        if xabs >= Self::from_uint(2u8) {
            /* |x| >= 2 */
            return xabs.bessel1_common(false, sign)
        }
        let half = Self::from(0.5);
        let mut z;
        if xabs >= Self::from((EXP_BASE as f64).powf(-13.0/127.0*Self::exp_bias().to_f64().unwrap()))
        {
            /* R0/S0 on [0,2] */
            const R0: [f64; 4] = [
                -6.25000000000000000000e-02, /* 0xBFB00000, 0x00000000 */
                1.40705666955189706048e-03, /* 0x3F570D9F, 0x98472C61 */
                -1.59955631084035597520e-05, /* 0xBEF0C5C6, 0xBA169668 */
                4.96727999609584448412e-08 /* 0x3E6AAAFA, 0x46CA0BD9 */
            ];
            const S0: [f64; 5] = [
                1.91537599538363460805e-02, /* 0x3F939D0B, 0x12637E53 */
                1.85946785588630915560e-04, /* 0x3F285F56, 0xB9CDF664 */
                1.17718464042623683263e-06, /* 0x3EB3BFF8, 0x333F8498 */
                5.04636257076217042715e-09, /* 0x3E35AC88, 0xC97DFF2C */
                1.23542274426137913908e-11 /* 0x3DAB2ACF, 0xCFB97ED8 */
            ];

            /* |x| >= 2**-127 */
            z = self*self;

            let mut b = false;
            let mut r = z*util::polynomial(&R0, z, SIGN_BIT, false);
            let mut s = Self::one() + z*util::polynomial(&S0, z, SIGN_BIT, false);
            if !SIGN_BIT && r.is_nan()
            {
                b = !b;
                r = z*util::polynomial(&R0.map(Neg::neg), z, SIGN_BIT, false);
            }
            if !SIGN_BIT && s.is_nan()
            {
                let one = Self::one();
                s = z*util::polynomial(&S0.map(Neg::neg), z, SIGN_BIT, false);
                if s > one
                {
                    b = !b;
                    s -= one
                }
                else
                {
                    s = one - s
                }
            }

            z = r/s;

            if b
            {
                if half < z
                {
                    return ((z - half)*self, true)
                }
                return ((half - z)*self, false)
            }
        }
        else
        {
            /* avoid underflow, raise inexact if x!=0 */
            z = self;
        }

        ((half + z)*self, false)
    }

    /// Bessel function of the first kind with α = 1.
    pub fn j1(self) -> Self
    {
        let (mut y, s) = self.j1_extra_sign();
        if s
        {
            y = -y
        }
        y
    }

    fn y1_extra_sign(self) -> (Self, bool)
    {
        /* y1(nan)=nan, y1(<0)=nan, y1(0)=-inf, y1(inf)=0 */
        if self.is_zero()
        {
            if !SIGN_BIT
            {
                return (Self::infinity(), true)
            }
            return (Self::neg_infinity(), false)
        }
        if self.is_sign_negative()
        {
            return (Self::snan(), false)
        }
        if self.is_infinite()
        {
            return (self.recip(), false)
        }

        if self >= Self::from_uint(2u8)
        {
            /* x >= 2 */
            return self.bessel1_common(true, false)
        }
        
        const TPI: f64 = 6.36619772367581382433e-01; /* 0x3FE45F30, 0x6DC9C883 */
        
        const U0: [f64; 5] = [
            -1.96057090646238940668e-01, /* 0xBFC91866, 0x143CBC8A */
            5.04438716639811282616e-02,  /* 0x3FA9D3C7, 0x76292CD1 */
            -1.91256895875763547298e-03, /* 0xBF5F55E5, 0x4844F50F */
            2.35252600561610495928e-05,  /* 0x3EF8AB03, 0x8FA6B88E */
            -9.19099158039878874504e-08, /* 0xBE78AC00, 0x569105B8 */
        ];
        const V0: [f64; 5] = [
            1.99167318236649903973e-02, /* 0x3F94650D, 0x3F4DA9F0 */
            2.02552581025135171496e-04, /* 0x3F2A8C89, 0x6C257764 */
            1.35608801097516229404e-06, /* 0x3EB6C05A, 0x894E8CA6 */
            6.22741452364621501295e-09, /* 0x3E3ABF1D, 0x5BA69A86 */
            1.66559246207992079114e-11, /* 0x3DB25039, 0xDACA772A */
        ];

        if self < Self::from((EXP_BASE as f64).powf(-3.125*EXP_SIZE as f64))
        {
            if !SIGN_BIT
            {
                return (Self::from(TPI)/self, true)
            }
            return (-Self::from(TPI)/self, false);
        }

        let one = Self::one();
        let z = self*self;

        let mut b = false;
        let mut u = util::polynomial(&U0, z, SIGN_BIT, false);
        let mut v = Self::one() + z*util::polynomial(&V0, z, SIGN_BIT, false);
        if !SIGN_BIT && u.is_nan()
        {
            b = !b;
            u = util::polynomial(&U0.map(Neg::neg), z, SIGN_BIT, false)
        }
        if !SIGN_BIT && v.is_nan()
        {
            v = z*util::polynomial(&V0.map(Neg::neg), z, SIGN_BIT, false);
            if v > one
            {
                b = !b;
                v -= one
            }
            else
            {
                v = one - v
            }
        }
        
        let uv = u/v;
        let xuv = self*uv;

        let (xj1, xj1_s) = self.j1_extra_sign();
        let xinv = self.recip();
        if b
        {
            if self < one
            {
                let xinvln = xinv.ln();
                let xj1mxinvln = xj1*xinvln;
                
                if xj1_s
                {
                    if xinv < xj1mxinvln
                    {
                        let o = Self::from(TPI)*(xj1mxinvln - xinv);
                        if xuv < o
                        {
                            return (o - xuv, false)
                        }
                        return (xuv - o, true)
                    }
                    return (Self::from(TPI)*(xinv - xj1mxinvln) + xuv, true)
                }

                return (Self::from(TPI)*(xj1mxinvln + xinv) + xuv, true)
            }
            let xln = self.ln();
            let xj1mxln = xj1*xln;

            if xj1_s
            {
                return (Self::from(TPI)*(xinv + xj1mxln) + xuv, true)
            }

            if xj1mxln < xinv
            {
                return (Self::from(TPI)*(xinv - xj1mxln) + xuv, true)
            }
            let o = Self::from(TPI)*(xj1mxln - xinv);
            if o < xuv
            {
                return (xuv - o, true)
            }
            return (o - xuv, false)
        }
        if !SIGN_BIT && self < one
        {
            if xj1_s
            {
                let xj1mxinvln = xj1*xinv.ln();
                if xinv < xj1mxinvln
                {
                    let o = Self::from(TPI)*(xj1mxinvln - xinv);
                    return (xuv + o, false)
                }
                let o = Self::from(TPI)*(xinv - xj1mxinvln);
                if xuv < o
                {
                    return (o - xuv, true)
                }
                return (xuv - o, false)
            }
            let o = Self::from(TPI)*(xj1*xinv.ln() + xinv);
            if xuv < o
            {
                return (o - xuv, true)
            }
            return (xuv - o, false)
        }
        let xln = self.ln();
        let xj1mxln = xj1*xln;
        if xj1_s
        {
            let o = Self::from(TPI)*(xinv + xj1mxln);
            if xuv < o
            {
                return (o - xuv, true)
            }
            return (xuv - o, false)
        }
        if xj1mxln < xinv
        {
            let o = Self::from(TPI)*(xinv - xj1mxln);
            if xuv < o
            {
                return (o - xuv, true)
            }
            return (xuv - o, false)
        }
        (xuv + Self::from(TPI)*(xj1mxln - xinv), false)
    }
    
    /// Bessel function of the second kind with α = 1.
    pub fn y1(self) -> Self
    {
        let (mut y, s) = self.y1_extra_sign();
        if s
        {
            y = -y
        }
        y
    }

    /// Bessel function of the first kind with α = `n`.
    pub fn jn(mut self, n: i32) -> Self
    {
        const INVSQRTPI: f64 = 5.64189583547756279280e-01; /* 0x3FE20DD7, 0x50429B6D */

        let mut sign = !self.sign_bit().is_zero();

        if self.is_nan()
        {
            /* nan */
            return self;
        }

        /* J(-n,x) = (-1)^n * J(n, x), J(n, -x) = (-1)^n * J(n, x)
        * Thus, J(-n,x) = J(n,-x)
        */
        /* nm1 = |n|-1 is used instead of |n| to handle n==INT_MIN */
        if n == 0
        {
            return self.j0();
        }
        let nm1;
        if n < 0
        {
            nm1 = -(n + 1);
            self = -self;
            sign = !sign;
        }
        else
        {
            nm1 = n - 1;
        }
        if nm1 == 0
        {
            return self.j1();
        }

        sign &= (n & 1) != 0; /* even n: 0, odd n: signbit(x) */
        self = self.abs();

        let one = Self::one();
        let two = Self::from_uint(2u8);

        let mut b_s;
        let mut b;
        if self.is_zero() || self.is_infinite()
        {
            /* if x is 0 or inf */
            b_s = false;
            b = Self::zero();
        }
        else if Self::from_int(nm1) < self
        {
            /* Safe to use J(n+1,x)=2n/x *J(n,x)-J(n-1,x) */
            if self >= Self::from(2.1359870359209100823950217061696e96)
            {
                /* x > 2**302 */
                /* (x >> n**2)
                *      Jn(x) = cos(x-(2n+1)*pi/4)*sqrt(2/x*pi)
                *      Yn(x) = sin(x-(2n+1)*pi/4)*sqrt(2/x*pi)
                *      Let s=sin(x), c=cos(x),
                *          xn=x-(2n+1)*pi/4, sqt2 = sqrt(2),then
                *
                *             n    sin(xn)*sqt2    cos(xn)*sqt2
                *          ----------------------------------
                *             0     s-c             c+s
                *             1    -s-c            -c+s
                *             2    -s+c            -c-s
                *             3     s+c             c-s
                */
                let (s, s_s) = self.sin_extra_sign();
                let (c, c_s) = self.cos_extra_sign();
                let (s1, s2) = match nm1 & 0b11 {
                    0b00 => (!c_s, s_s),
                    0b01 => (!c_s, !s_s),
                    0b10 => (c_s, !s_s),
                    0b11 => (c_s, s_s),
                    _ => unreachable!()
                };
                let temp;
                (temp, b_s) = util::add_extra_sign((c, s1), (s, s2), SIGN_BIT);
                b = Self::from(INVSQRTPI)*temp/self.sqrt();
            }
            else
            {
                let (mut a, mut a_s) = self.j0_extra_sign();
                (b, b_s) = self.j1_extra_sign();
                let mut i = 0;
                while i < nm1 {
                    i += 1;
                    let (temp, temp_s) = (b, b_s);
                    b *= two*Self::from_int(i)/self;
                    (b, b_s) = util::add_extra_sign((b, b_s), (a, !a_s), SIGN_BIT); /* avoid underflow */
                    a_s = temp_s;
                    a = temp;
                }
            }
        }
        else if self < Self::from(0.00000000186264514923095703125)
        {
            /* x < 2**-29 */
            /* x is tiny, return the first Taylor expansion of J(n,x)
            * J(n,x) = 1/n!*(x/2)^n  - ...
            */
            if nm1 > 32
            {
                /* underflow */
                b_s = false;
                b = Self::zero();
            }
            else
            {
                let temp = self*Self::from(0.5);
                b_s = false;
                b = temp;
                let mut a = one;
                let mut i = 2;
                while i <= nm1 + 1
                {
                    a *= Self::from_int(i); /* a = n! */
                    b *= temp; /* b = (x/2)^n */
                    i += 1;
                }
                b /= a;
            }
        }
        else
        {
            /* use backward recurrence */
            /*                      x      x^2      x^2
            *  J(n,x)/J(n-1,x) =  ----   ------   ------   .....
            *                      2n  - 2(n+1) - 2(n+2)
            *
            *                      1      1        1
            *  (for large x)   =  ----  ------   ------   .....
            *                      2n   2(n+1)   2(n+2)
            *                      -- - ------ - ------ -
            *                       x     x         x
            *
            * Let w = 2n/x and h=2/x, then the above quotient
            * is equal to the continued fraction:
            *                  1
            *      = -----------------------
            *                     1
            *         w - -----------------
            *                        1
            *              w+h - ---------
            *                     w+2h - ...
            *
            * To determine how many terms needed, let
            * Q(0) = w, Q(1) = w(w+h) - 1,
            * Q(k) = (w+k*h)*Q(k-1) - Q(k-2),
            * When Q(k) > 1e4      good for single
            * When Q(k) > 1e9      good for double
            * When Q(k) > 1e17     good for quadruple
            */
            /* determine k */

            let nf = Self::from_int(nm1) + one;
            let w = two*nf/self;
            let h = two/self;
            let mut z = w + h;
            let wz = w*z;

            let mut q0_s = false;
            let mut q0 = w;
            let mut q1_s = false;
            let mut q1 = if SIGN_BIT || wz >= one
            {
                wz - one
            }
            else
            {
                q1_s = !q1_s;
                one - wz
            };

            let mut k = 1;
            let q1_min = Self::from((EXP_BASE as f64).powf(13.287712379549449391481277717958/8.0*EXP_SIZE as f64));
            while q1 < q1_min
            {
                k += 1;
                z += h;
                let (tmp, tmp_s) = util::add_extra_sign((z*q1, q1_s), (q0, !q0_s), SIGN_BIT);
                q0 = q1;
                q0_s = q1_s;
                q1 = tmp;
                q1_s = tmp_s;
            }
            let mut t = Self::zero();
            let mut i = k;
            while i >= 0
            {
                t = (two*(Self::from_int(i) + nf)/self - t).recip();
                assert!(!t.is_nan());
                i -= 1;
            }
            let mut a_s = false;
            let mut a = t;
            b_s = false;
            b = one;
            /*  estimate log((2/x)^n*n!) = n*log(2/x)+n*ln(n)
            *  Hence, if n*(log(2n/x)) > ...
            *  single 8.8722839355e+01
            *  double 7.09782712893383973096e+02
            *  long double 1.1356523406294143949491931077970765006170e+04
            *  then recurrent value may overflow and the result is
            *  likely underflow to zero
            */
            let wabs = w.abs();
            let tmp = nf*wabs.ln();
            if tmp < Self::from(0.69860503429133858267716535433071*Self::exp_bias().to_f64().unwrap()) || wabs < one
            {
                i = nm1;
                while i > 0
                {
                    let (temp, temp_s) = (b, b_s);
                    b *= (two*Self::from_int(i))/self;
                    (b, b_s) = util::add_extra_sign((b, b_s), (a, !a_s), SIGN_BIT);
                    a = temp;
                    a_s = temp_s;
                    i -= 1;
                }
            }
            else
            {
                let x1p500 = Self::from((EXP_BASE as f64).powf(0.47244094488188976377952755905512*Self::exp_bias().to_f64().unwrap())); // 0x1p500 == 2^500

                i = nm1;
                while i > 0
                {
                    let (temp, temp_s) = (b, b_s);
                    b *= (two*Self::from_int(i))/self;
                    (b, b_s) = util::add_extra_sign((b, b_s), (a, !a_s), SIGN_BIT);
                    a = temp;
                    a_s = temp_s;

                    /* scale b to avoid spurious overflow */
                    if b > x1p500
                    {
                        a /= b;
                        t /= b;
                        b = one;
                    }
                    i -= 1;
                }
            }
            let (z, z_s) = self.j0_extra_sign();
            let (w, w_s) = self.j1_extra_sign();
            if z.abs() >= w.abs()
            {
                b_s ^= z_s;
                b = t * z / b;
            }
            else
            {
                b_s = w_s^a_s;
                b = t * w / a;
            }
        }

        if sign^b_s
        {
            b = -b
        }
        
        b
    }

    /// Bessel function of the second kind with α = `n`.
    pub fn yn(self, n: i32) -> Self
    {
        const INVSQRTPI: f64 = 5.64189583547756279280e-01; /* 0x3FE20DD7, 0x50429B6D */

        if self.is_nan()
        {
            /* nan */
            return self;
        }
        let zero = Self::zero();
        if self < zero
        {
            /* x < 0 */
            return Self::snan();
        }
        if self.is_infinite()
        {
            return zero;
        }
    
        if n == 0
        {
            return self.y0();
        }
        let nm1;
        let sign;
        if n < 0
        {
            nm1 = -(n + 1);
            sign = (n & 1) != 0;
        }
        else
        {
            nm1 = n - 1;
            sign = false;
        }
        if nm1 == 0
        {
            let (mut y1, y1_s) = self.y1_extra_sign();

            if sign^y1_s
            {
                y1 = -y1
            }

            return y1
        }
    
        let two = Self::from(2u8);

        let mut b_s;
        let mut b;
        if self > Self::from(8.1481439053379443450737827536375e90)
        {
            /* x > 2**302 */
            /* (x >> n**2)
             *      Jn(x) = cos(x-(2n+1)*pi/4)*sqrt(2/x*pi)
             *      Yn(x) = sin(x-(2n+1)*pi/4)*sqrt(2/x*pi)
             *      Let s=sin(x), c=cos(x),
             *          xn=x-(2n+1)*pi/4, sqt2 = sqrt(2),then
             *
             *             n    sin(xn)*sqt2    cos(xn)*sqt2
             *          ----------------------------------
             *             0     s-c             c+s
             *             1    -s-c            -c+s
             *             2    -s+c            -c-s
             *             3     s+c             c-s
             */
            let (s, s_s) = self.sin_extra_sign();
            let (c, c_s) = self.cos_extra_sign();
            let (s1, s2) = match nm1 & 0b11 {
                0b00 => (!s_s, !c_s),
                0b01 => (!s_s, c_s),
                0b10 => (s_s, c_s),
                0b11 => (s_s, !c_s),
                _ => unreachable!()
            };
            let temp;
            (temp, b_s) = util::add_extra_sign((s, s1), (c, s2), SIGN_BIT);
            b = Self::from(INVSQRTPI)*temp/self.sqrt();
        }
        else
        {
            let (mut a, mut a_s) = self.y0_extra_sign();
            (b, b_s) = self.y1_extra_sign();
            /* quit if b is -inf */
            let mut i = 0;
            while i < nm1 && b.is_finite()
            {
                i += 1;
                let (temp, temp_s) = (b, b_s);
                b *= two*Self::from(i)/self;
                (b, b_s) = util::add_extra_sign((b, b_s), (a, !a_s), SIGN_BIT);
                (a, a_s) = (temp, temp_s);
            }
        }
    
        if sign^b_s
        {
            b = -b
        }

        b
    }
}

impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:,
    U: ConstZero
{
    /// The additive identity element of `Self`, `0`.
    pub const ZERO: Self = Self::from_bits(U::ZERO);
}

#[cfg(test)]
mod test
{
    #![allow(unused)]

    use crate::{g_711::FpG711, ieee754::{FpDouble, FpHalf}, tests::F, Fp};

    #[test]
    fn test_gamma()
    {
        crate::tests::test_op1("ln_gamma", |x| f32::ln_gamma(x).0, |x| Fp::ln_gamma(x).0, None, Some(-(crate::tests::F::SIGN_SIZE as f32)*4.5..20.0));
        crate::tests::test_op1("gamma", f32::gamma, Fp::gamma, None, Some(-(crate::tests::F::SIGN_SIZE as f32)*4.5..7.0))
    }
    
    #[test]
    fn test_j0()
    {
        crate::tests::test_op1("j0", libm::j0f, Fp::j0, Some(0.1), Some(-(crate::tests::F::SIGN_SIZE as f32)*20.0..20.0))
    }
    
    #[test]
    fn test_y0()
    {
        //crate::tests::test_op1("y0", libm::y0f, |x| crate::tests::F::from(x.y0_extra_sign().1 as u8), Some(0.1), Some(0.01..20.0))
        //crate::tests::test_op1("y0", libm::y0f, |x| x.y0_extra_sign().0, Some(0.1), Some(0.01..20.0))
        crate::tests::test_op1("y0", libm::y0f, Fp::y0, Some(0.1), Some(0.01..20.0))
    }
    
    #[test]
    fn test_j1()
    {
        crate::tests::test_op1("j1", libm::j1f, Fp::j1, Some(0.1), Some(-(crate::tests::F::SIGN_SIZE as f32)*20.0..20.0))
    }
    
    #[test]
    fn test_y1()
    {
        crate::tests::test_op1("y1", libm::y1f, Fp::y1, Some(0.1), Some(0.1..20.0))
    }
    
    #[test]
    fn test_j2()
    {
        crate::tests::test_op1("j2", |x| libm::jnf(2, x), |x| x.jn(2), Some(0.1), Some(-(crate::tests::F::SIGN_SIZE as f32)*20.0..20.0))
    }
    
    #[test]
    fn test_y2()
    {
        crate::tests::test_op1("y2", |x| libm::ynf(2, x), |x| x.yn(2), Some(0.1), Some(1.0..20.0))
    }
    
    #[test]
    fn test_j3()
    {
        crate::tests::test_op1("j3", |x| libm::jnf(3, x), |x| x.jn(3), Some(0.1), Some(-(crate::tests::F::SIGN_SIZE as f32)*20.0..20.0))
    }
    
    #[test]
    fn test_y3()
    {
        crate::tests::test_op1("y3", |x| libm::ynf(3, x), |x| x.yn(3), Some(0.1), Some(1.0..20.0))
    }
    
    #[test]
    fn test_erf()
    {
        crate::tests::test_op1("erf", libm::erff, Fp::erf, Some(0.1), Some(-(crate::tests::F::SIGN_SIZE as f32)*5.0..5.0))
    }

    #[test]
    fn test_erfc_once()
    {
        let x = F::from(0.9);

        let y = x.erfc();

        println!("{y}")
    }
    
    #[test]
    fn test_erfc()
    {
        crate::tests::test_op1("erfc", libm::erfcf, Fp::erfc, Some(0.1), Some(-(crate::tests::F::SIGN_SIZE as f32)*5.0..5.0))
    }

    #[test]
    fn test_next_up_down()
    {
        type F = Fp<u8, true, 3, 1, 3, {usize::MAX}>;

        let mut x = F::neg_infinity();

        loop
        {
            let y = x.next_up();
            if !y.is_zero()
            {
                if x == y
                {
                    let y = x.next_up();
                }
                assert_ne!(x, y);
            }
            assert_eq!(-(-x).next_down(), y);
            if !(x == y.next_down())
            {
                let yy = x.next_up();
                let xx = y.next_down();
                println!("{x} ^ {yy} v {xx}")
            }
            assert_eq!(x, y.next_down());
            x = y;
            if !x.is_finite()
            {
                break
            }
        }

        loop
        {
            let y = x.next_down();
            if !y.is_zero()
            {
                assert_ne!(x, y);
            }
            assert_eq!(-(-x).next_up(), y);
            assert_eq!(x, y.next_up());
            x = y;
            if !x.is_finite()
            {
                break
            }
        }
    }
}