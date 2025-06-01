#![allow(unused)]

use core::ops::{Add, AddAssign, Div, MulAssign, Neg, Rem, Shl, Shr};

use num_traits::{CheckedShl, CheckedAdd, CheckedSub, Float, NumCast, One, Zero};

use crate::{ieee754::{FpDouble, FpHalf, FpQuadruple, FpSingle}, AnyInt, Fp, Int, UInt};

pub(crate) fn do_nothing<T>(x: T) -> T
{
    x
}

#[derive(Clone, Copy)]
pub(crate) enum Saturation
{
    Overflow,
    Underflow
}

impl Neg for Saturation
{
    type Output = Self;

    fn neg(self) -> Self::Output
    {
        match self
        {
            Self::Overflow => Self::Underflow,
            Self::Underflow => Self::Overflow
        }
    }
}

pub(crate) trait Conversion//: Float
{
    type Bits: UInt;

    /// Wether or not the number has a sign bit. If not, it can only be positive.
    const SIGN_BIT: bool;
    /// The size of the exponent part in bits
    const EXP_SIZE: usize;
    /// The size of the integer part in bits
    const INT_SIZE: usize;
    /// The size of the fractional part in bits
    const FRAC_SIZE: usize;
    /// The base for the exponent
    const EXP_BASE: usize;
}

impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Conversion for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:
{
    type Bits = U;
    
    const SIGN_BIT: bool = SIGN_BIT;
    const EXP_SIZE: usize = EXP_SIZE;
    const INT_SIZE: usize = INT_SIZE;
    const FRAC_SIZE: usize = FRAC_SIZE;
    const EXP_BASE: usize = EXP_BASE;
}

macro_rules! impl_float {
    ($($f:ty: $fp:ty),*) => {
        $(
            impl Conversion for $f
            {
                type Bits = <$fp as Conversion>::Bits;
            
                const SIGN_BIT: bool = <$fp as Conversion>::SIGN_BIT;
                const EXP_SIZE: usize = <$fp as Conversion>::EXP_SIZE;
                const INT_SIZE: usize = <$fp as Conversion>::INT_SIZE;
                const FRAC_SIZE: usize = <$fp as Conversion>::FRAC_SIZE;
                const EXP_BASE: usize = <$fp as Conversion>::EXP_BASE;
            }
        )*
    };
}
impl_float!(f16: FpHalf, f32: FpSingle, f64: FpDouble, f128: FpQuadruple);

macro_rules! impl_int {
    ($($i:ty: $u:ty),*) => {
        $(
            impl Conversion for $u
            {
                type Bits = $u;
            
                const SIGN_BIT: bool = false;
                const EXP_SIZE: usize = 0;
                const INT_SIZE: usize = bitsize_of::<$u>();
                const FRAC_SIZE: usize = 0;
                const EXP_BASE: usize = 2;
            }
            impl Conversion for $i
            {
                type Bits = $u;
            
                const SIGN_BIT: bool = true;
                const EXP_SIZE: usize = 0;
                const INT_SIZE: usize = bitsize_of::<$u>() - 1;
                const FRAC_SIZE: usize = 0;
                const EXP_BASE: usize = 2;
            }
        )*
    };
}
impl_int!(i8: u8, i16: u16, i32: u32, i64: u64, i128: u128);

pub const fn is_float_conversion_mutually_lossless<A, B>() -> bool
where
    A: Conversion,
    B: Conversion
{
    is_float_conversion_lossless::<A, B>() && is_float_conversion_lossless::<B, A>()
}

pub const fn is_float_conversion_lossless<From, To>() -> bool
where
    From: Conversion,
    To: Conversion
{
    if From::SIGN_BIT && !To::SIGN_BIT
        || From::INT_SIZE > To::INT_SIZE
        || From::FRAC_SIZE > To::FRAC_SIZE
    {
        return false
    }
    if To::EXP_BASE == 0
    {
        To::EXP_SIZE == 0
            && From::EXP_SIZE == 0
    }
    else
    {
        if From::EXP_BASE < To::EXP_BASE
            || From::EXP_BASE % To::EXP_BASE != 0
        {
            return false
        }
        let mut b = From::EXP_BASE/To::EXP_BASE;
        let mut a = To::EXP_SIZE;
        while b != 0 && a % 2 == 0
        {
            a >>= 1;
            b -= 1
        }
        b <= u32::MAX as usize || match From::EXP_SIZE.checked_shl(b as u32)
        {
            Some(b) => b <= a,
            None => false
        }
    }
}

trait AnyIntSpec: AnyInt
{
    const IS_SIGNED: bool;
}
impl<T> AnyIntSpec for T
where
    T: AnyInt
{
    default const IS_SIGNED: bool = true;
}
impl<T> AnyIntSpec for T
where
    T: UInt
{
    const IS_SIGNED: bool = false;
}

pub const fn is_signed<T: AnyInt>() -> bool
{
    <T as AnyIntSpec>::IS_SIGNED
}

pub fn utilized_size_of<T: AnyInt>(x: &T) -> usize
{
    bitsize_of::<T>() - x.leading_zeros() as usize + is_signed::<T>() as usize
}

pub const fn icbrt(a: usize) -> usize
{
    if bitsize_of::<usize>() <= 32
    {
        // Implementation based on Hacker's Delight `icbrt2`
        let mut x = a;
        let mut y2 = 0;
        let mut y = 0;
        let smax = bitsize_of::<usize>() / 3;
        let mut s = smax + 1;
        while s > 0
        {
            s -= 1;
            let s = s * 3;
            y2 *= 4;
            y *= 2;
            let b = 3 * (y2 + y) + 1;
            if x >> s >= b {
                x -= b << s;
                y2 += 2 * y + 1;
                y += 1;
            }
        }
        return y;
    }

    if a < 8
    {
        return (a > 0) as usize;
    }

    #[inline]
    const fn guess(x: usize) -> usize
    {
        1 << ((x.ilog2() + 2) / 3)
    }

    #[inline]
    const fn fixpoint(a: usize, mut x: usize) -> usize
    {
        let mut xn = (a / (x * x) + x * 2) / 3;
        while x < xn {
            x = xn;
            xn = (a / (x * x) + x * 2) / 3;
        }
        while x > xn {
            x = xn;
            xn = (a / (x * x) + x * 2) / 3;
        }
        x
    }

    fixpoint(a, guess(a))
}

pub const fn bitsize_of<T>() -> usize
{
    core::mem::size_of::<T>()*8
}

pub fn complementary_add_sub_assign<T>(y_add: &mut T, y_sub: &mut T, mut x: T) -> Result<(), ()>
where
    T: CheckedSub + CheckedAdd + Copy
{
    if let Some(diff) = y_sub.checked_sub(&x)
    {
        *y_sub = diff;
        return Ok(())
    }
    x = x - *y_sub;
    if let Some(sum) = y_add.checked_add(&x)
    {
        *y_add = sum;
        return Ok(())
    }
    Err(())
}

#[inline]
pub fn polynomial<T, F, const N: usize>(p: &[T; N], mut z: F, may_be_neg: bool, z_neg: bool) -> F
where
    T: Float,
    F: Float + From<T> + AddAssign + MulAssign
{
    if z_neg && may_be_neg
    {
        z = -z;
    }
    let mut y = F::zero();
    let mut zn = F::one();
    for (n, p) in p.iter()
        .enumerate()
    {
        if may_be_neg
        {
            y += <F as From<_>>::from(*p)*zn;
        }
        else if p.is_sign_positive() ^ (z_neg && (n % 2 != 0))
        {
            y += <F as From<_>>::from(p.abs())*zn;
        }
        zn *= z;
    }
    if !may_be_neg
    {
        let mut zn = F::one();
        for (n, p) in p.iter()
            .enumerate()
        {
            if p.is_sign_negative() ^ (z_neg && (n % 2 != 0))
            {
                y = y - <F as From<_>>::from(p.abs())*zn;
            }
            zn *= z;
        }
    }
    y
}

pub fn chebychev_approximation<T, const N: usize>(coeffs: [T; N]) -> [T; N]
where
    T: Float + AddAssign
{
    let t = (0..N).map(
        |n| chebyshev_polynomial::<T, N>(1, n).unwrap()
    );
    t.zip(coeffs)
        .map(|(t, c)| t.map(|tn| c*tn))
        .reduce(|a, b| unsafe {
            a.into_iter()
                .zip(b)
                .map(|(a, b)| a + b)
                .next_chunk()
                .unwrap_unchecked()
        }).unwrap_or_else(|| [T::zero(); N])
}

pub fn chebyshev_polynomial<T, const N: usize>(kind: usize, order: usize) -> Option<[T; N]>
where
    T: Float + AddAssign
{
    if order > N
    {
        return None
    }

    let two = T::one() + T::one();
    let mut t_prev = [T::zero(); N];
    t_prev[0] = T::one();
    if order == 0
    {
        return Some(t_prev)
    }
    
    let mut kind_c = T::zero();
    let mut k = 0;
    while k < kind
    {
        kind_c += T::one();
        k += 1;
    }

    let mut t: [T; N] = core::array::from_fn(|i| if i == 1 {kind_c} else {T::zero()});

    let mut k = 1;
    while k < order
    {
        let mut t_prev_iter = t_prev.into_iter();
        let mut t_iter = t.into_iter();
        let mut first = true;
        
        let t_next = core::array::from_fn(|_| if first
            {
                first = false;
                -t_prev_iter.next().unwrap()
            }
            else
            {
                two * t_iter.next().unwrap() - t_prev_iter.next().unwrap()
            }
        );

        t_prev = t;
        t = t_next;
        k += 1;
    }

    Some(t)
}

/// Returns (2^`exponent`).ilog(`base`)
pub const fn exp2_ilog(exponent: usize, base: usize) -> usize
{
    if base == 0
    {
        return 0
    }
    if base == 2
    {
        return exponent
    }
    if base.is_power_of_two()
    {
        return exponent/base.ilog2() as usize
    }
    let mut d: usize = 1;
    let mut o = 0;

    let mut n = exponent;
    while n > 0
    {
        if d.leading_zeros() == 0
        {
            o += 1;
            d /= base;
        }
        d <<= 1;
        n -= 1
    }

    let y = o + d.ilog(base) as usize;

    #[cfg(debug_assertions)]
    if let Some(s) = 1usize.checked_shl(exponent as u32)
    {
        debug_assert!(y == s.ilog(base) as usize)
    }

    y
}

/// Returns (`base`^`exponent`).ilog2()
pub const fn pow_ilog2(mut exponent: usize, base: usize) -> usize
{
    if base == 0 || base == 1
    {
        return 0
    }
    if base.is_power_of_two()
    {
        return exponent*base.ilog2() as usize
    }
    let p = base.ilog2();
    let mut d: usize = 1;
    let mut o = 0;

    let mut n = exponent;
    while n > 0
    {
        if d.leading_zeros() <= p
        {
            o += 1;
            d >>= 1;
        }
        d *= base;
        n -= 1
    }

    let y = o + d.ilog2() as usize;

    #[cfg(debug_assertions)]
    if let Some(s) = base.checked_pow(exponent as u32)
    {
        debug_assert!(y == s.ilog2() as usize)
    }

    y
}

pub const fn is_power_of(mut value: usize, base: usize) -> bool
{
    if value == 0
    {
        return false
    }

    loop
    {
        if value == 1
        {
            return true
        }
        if !value.is_multiple_of(base)
        {
            return false
        }
        value /= base
    }
}

trait RoundingDivPowSpec<I: UInt>: Div<Output = Self> + Rem<Output = Self> + Shl<u32, Output = Self> + Add<Output = Self> + PartialOrd + One + Copy
{
    fn _rounding_div_pow(x: Self, rhs: Self, pow: I) -> Self;
}
impl<T, I: UInt> RoundingDivPowSpec<I> for T
where
    Self: Div<Output = Self> + Rem<Output = Self> + Shl<u32, Output = Self> + Add<Output = Self> + PartialOrd + One + Copy
{
    default fn _rounding_div_pow(x: Self, rhs: Self, pow: I) -> Self
    {
        rounding_div_pow_naive(x, rhs, pow)
    }
}
impl<T, I: UInt> RoundingDivPowSpec<I> for T
where
    Self: Div<Output = Self> + Rem<Output = Self> + Shl<u32, Output = Self> + Shr<u32, Output = Self> + Add<Output = Self> + PartialOrd + One + Copy
    + NumCast + CheckedShl
{
    fn _rounding_div_pow(x: Self, rhs: Self, pow: I) -> Self
    {
        if T::from(2usize) == Some(rhs) && let Some(pow) = pow.to_u32() && let Some(rhs) = rhs.checked_shl(pow)
        {
            return rounding_div_2_pow(x, pow)
        }
        rounding_div_pow_naive(x, rhs, pow)
    }
}

fn rounding_div_pow_naive<T, I: UInt>(mut x: T, rhs: T, mut pow: I) -> T
where
    T: Div<T, Output = T> + Rem<T, Output = T> + Shl<u32, Output = T> + Add<T, Output = T> + PartialOrd + One + Copy
{
    while pow > I::zero()
    {
        x = rounding_div(x, rhs);
        pow = pow - I::one();
    }
    x
}

pub fn rounding_div_pow<T, I: UInt>(mut x: T, rhs: T, mut pow: I) -> T
where
    T: Div<T, Output = T> + Rem<T, Output = T> + Shl<u32, Output = T> + Add<T, Output = T> + PartialOrd + One + Copy
{
    T::_rounding_div_pow(x, rhs, pow)
}

pub fn rounding_div<T>(x: T, rhs: T) -> T
where
    T: Div<T, Output = T> + Rem<T, Output = T> + Shl<u32, Output = T> + Add<T, Output = T> + PartialOrd + One + Copy
{
    if (x % rhs) << 1 > rhs
    {
        x/rhs + T::one()
    }
    else
    {
        x/rhs
    }
}

pub fn ceil_div<T>(x: T, rhs: T) -> T
where
    T: Div<T, Output = T> + Rem<T, Output = T> + Shl<u32, Output = T> + Add<T, Output = T> + PartialOrd + One + Copy + Zero
{
    if (x % rhs) > T::zero()
    {
        x/rhs + T::one()
    }
    else
    {
        x/rhs
    }
}

pub fn rounding_div_2_pow<T, I: UInt>(x: T, pow: I) -> T
where
    T: Shr<u32, Output = T> + Rem<T, Output = T> + Shl<u32, Output = T> + Add<T, Output = T> + NumCast + PartialOrd + One + Copy + Shr<I, Output = T> + Shl<I, Output = T>
{
    let two = T::one() << pow;
    if (x % two) << 1 > two
    {
        (x >> pow) + T::one()
    }
    else
    {
        x >> pow
    }
}

pub fn rounding_div_2<T>(x: T) -> T
where
    T: Shr<u32, Output = T> + Rem<T, Output = T> + Shl<u32, Output = T> + Add<T, Output = T> + NumCast + PartialOrd + One + Copy
{
    let two = T::from(2).unwrap();
    if (x % two) << 1 > two
    {
        (x >> 1) + T::one()
    }
    else
    {
        x >> 1
    }
}

pub fn powi<P, I: Int>(x: P, mut n: I) -> P
where
    P: One + Div<P, Output = P> + Copy
{
    let mut r = P::one();
    
    let mut x = if n < I::zero()
    {
        n = -n;
        P::one() / x
    }
    else
    {
        x
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

pub fn powu<P, I: UInt>(mut x: P, mut n: I) -> P
where
    P: One + Copy
{
    let mut r = P::one();
    
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
