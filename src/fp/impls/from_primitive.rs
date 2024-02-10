use num_traits::FromPrimitive;

use crate::{bitsize_of, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const INT_BIT: bool, const FRAC_SIZE: usize> FromPrimitive for Fp<U, EXP_SIZE, INT_BIT, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - INT_BIT as usize - FRAC_SIZE - 1]:,
    [(); bitsize_of::<U>() - EXP_SIZE - false as usize - FRAC_SIZE - 1]:
{
    fn from_isize(n: isize) -> Option<Self>
    {
        Some(Self::from_int(n))
    }

    fn from_i8(n: i8) -> Option<Self>
    {
        Some(Self::from_int(n))
    }

    fn from_i16(n: i16) -> Option<Self>
    {
        Some(Self::from_int(n))
    }

    fn from_i32(n: i32) -> Option<Self>
    {
        Some(Self::from_int(n))
    }

    fn from_i64(n: i64) -> Option<Self>
    {
        Some(Self::from_int(n))
    }

    fn from_i128(n: i128) -> Option<Self>
    {
        Some(Self::from_int(n))
    }

    fn from_usize(n: usize) -> Option<Self>
    {
        Some(Self::from_uint(n))
    }

    fn from_u8(n: u8) -> Option<Self>
    {
        Some(Self::from_uint(n))
    }

    fn from_u16(n: u16) -> Option<Self>
    {
        Some(Self::from_uint(n))
    }

    fn from_u32(n: u32) -> Option<Self>
    {
        Some(Self::from_uint(n))
    }

    fn from_u64(n: u64) -> Option<Self>
    {
        Some(Self::from_uint(n))
    }

    fn from_u128(n: u128) -> Option<Self>
    {
        Some(Self::from_uint(n))
    }

    fn from_f32(n: f32) -> Option<Self>
    {
        Some(<Self as From<_>>::from(n))
    }

    fn from_f64(n: f64) -> Option<Self>
    {
        Some(<Self as From<_>>::from(n))
    }
}