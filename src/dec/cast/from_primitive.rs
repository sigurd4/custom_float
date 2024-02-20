use num_traits::FromPrimitive;

use crate::{util, Dec, UInt};

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> FromPrimitive for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:
{
    #[inline]
    fn from_isize(n: isize) -> Option<Self>
    {
        Some(Self::from_int(n))
    }

    #[inline]
    fn from_i8(n: i8) -> Option<Self>
    {
        Some(Self::from_int(n))
    }

    #[inline]
    fn from_i16(n: i16) -> Option<Self>
    {
        Some(Self::from_int(n))
    }

    #[inline]
    fn from_i32(n: i32) -> Option<Self>
    {
        Some(Self::from_int(n))
    }

    #[inline]
    fn from_i64(n: i64) -> Option<Self>
    {
        Some(Self::from_int(n))
    }

    #[inline]
    fn from_i128(n: i128) -> Option<Self>
    {
        Some(Self::from_int(n))
    }

    #[inline]
    fn from_usize(n: usize) -> Option<Self>
    {
        Some(Self::from_uint(n))
    }

    #[inline]
    fn from_u8(n: u8) -> Option<Self>
    {
        Some(Self::from_uint(n))
    }

    #[inline]
    fn from_u16(n: u16) -> Option<Self>
    {
        Some(Self::from_uint(n))
    }

    #[inline]
    fn from_u32(n: u32) -> Option<Self>
    {
        Some(Self::from_uint(n))
    }

    #[inline]
    fn from_u64(n: u64) -> Option<Self>
    {
        Some(Self::from_uint(n))
    }

    #[inline]
    fn from_u128(n: u128) -> Option<Self>
    {
        Some(Self::from_uint(n))
    }

    #[inline]
    fn from_f32(n: f32) -> Option<Self>
    {
        Some(<Self as From<_>>::from(n))
    }

    #[inline]
    fn from_f64(n: f64) -> Option<Self>
    {
        Some(<Self as From<_>>::from(n))
    }
}