#![allow(unused)]

use core::ops::{Add, Div, Rem, Shl, Shr};
use std::ops::{AddAssign, MulAssign};

use array_math::ArrayMath;
use num_traits::{NumCast, One, Zero, Float};

use crate::{Int, UInt};

pub const fn bitsize_of<T>() -> usize
{
    core::mem::size_of::<T>()*8
}

#[inline]
pub fn add_extra_sign<T>(a: (T, bool), b: (T, bool), may_be_neg: bool) -> (T, bool)
where
    T: Float
{
    let (a, a_s) = a;
    let (b, b_s) = b;
    let mut s = false;
    let y = match (a_s, b_s)
    {
        (false, false) => a + b,
        (true, false) => if may_be_neg || a <= b
        {
            b - a
        }
        else
        {
            s = !s;
            a - b
        },
        (false, true) => if may_be_neg || b <= a
        {
            a - b
        }
        else
        {
            s = !s;
            b - a
        },
        (true, true) => if may_be_neg
        {
            -a - b
        }
        else
        {
            s = !s;
            a + b
        }
    };

    (y, s)
}

#[inline]
pub fn polynomial<T, F, const N: usize>(p: &[T; N], z: F, may_be_neg: bool, z_neg: bool) -> F
where
    T: Float,
    F: Float + From<T> + AddAssign + MulAssign
{
    if may_be_neg
    {
        return p.map(<F as From<_>>::from).polynomial(if z_neg {-z} else {z})
    }
    let mut y = F::zero();
    let mut zn = F::one();
    for n in 0..N
    {
        if may_be_neg || p[n].is_sign_positive() ^ (z_neg && (n % 2 != 0))
        {
            y = y + <F as From<_>>::from(p[n].abs())*zn;
        }
        zn = zn*z;
    }
    if !may_be_neg
    {
        let mut zn = F::one();
        for n in 0..N
        {
            if p[n].is_sign_negative() ^ (z_neg && (n % 2 != 0))
            {
                y = y - <F as From<_>>::from(p[n].abs())*zn;
            }
            zn = zn*z;
        }
    }
    y
}

/// Returns (2^`digits_bin`).ilog(base)
pub fn count_digits_in_base(digits_bin: usize, base: usize) -> usize
{
    let mut d: usize = 1;
    let mut o = 0;

    for _ in 0..digits_bin
    {
        if d.leading_zeros() == 0
        {
            o += 1;
            d /= base;
        }
        d <<= 1
    }

    let y = o + d.ilog(base) as usize;

    if let Some(s) = 1usize.checked_shl(digits_bin as u32)
    {
        assert_eq!(y, s.ilog(base) as usize)
    }

    y
}

pub fn rounding_div_pow<T, I: UInt>(mut x: T, rhs: T, mut pow: I) -> T
where
    T: Div<T, Output = T> + Rem<T, Output = T> + Shl<usize, Output = T> + Add<T, Output = T> + PartialOrd + One + Copy
{
    while pow > I::zero()
    {
        x = rounding_div(x, rhs);
        pow = pow - I::one();
    }
    x
}

pub fn rounding_div<T>(x: T, rhs: T) -> T
where
    T: Div<T, Output = T> + Rem<T, Output = T> + Shl<usize, Output = T> + Add<T, Output = T> + PartialOrd + One + Copy
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
    T: Div<T, Output = T> + Rem<T, Output = T> + Shl<usize, Output = T> + Add<T, Output = T> + PartialOrd + One + Copy + Zero
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

pub fn rounding_div_2<T>(x: T) -> T
where
    T: Shr<usize, Output = T> + Rem<T, Output = T> + Shl<usize, Output = T> + Add<T, Output = T> + NumCast + PartialOrd + One + Copy
{
    let two = T::from(2).unwrap();
    if (x % two) << 1 > two
    {
        (x >> 1usize) + T::one()
    }
    else
    {
        x >> 1usize
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