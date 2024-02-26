use num_traits::ToPrimitive;

use crate::{util, Fp, UInt};

impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> ToPrimitive for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - 0 - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:
{
    #[inline]
    fn to_i8(&self) -> Option<i8>
    {
        self.to_int()
    }
    #[inline]
    fn to_u8(&self) -> Option<u8>
    {
        self.to_uint()
    }

    #[inline]
    fn to_i16(&self) -> Option<i16>
    {
        self.to_int()
    }
    #[inline]
    fn to_u16(&self) -> Option<u16>
    {
        self.to_uint()
    }

    #[inline]
    fn to_i32(&self) -> Option<i32>
    {
        self.to_int()
    }
    #[inline]
    fn to_u32(&self) -> Option<u32>
    {
        self.to_uint()
    }

    #[inline]
    fn to_isize(&self) -> Option<isize>
    {
        self.to_int()
    }
    #[inline]
    fn to_usize(&self) -> Option<usize>
    {
        self.to_uint()
    }

    #[inline]
    fn to_i64(&self) -> Option<i64>
    {
        self.to_int()
    }
    #[inline]
    fn to_u64(&self) -> Option<u64>
    {
        self.to_uint()
    }

    #[inline]
    fn to_i128(&self) -> Option<i128>
    {
        self.to_int()
    }
    #[inline]
    fn to_u128(&self) -> Option<u128>
    {
        self.to_uint()
    }

    #[inline]
    fn to_f32(&self) -> Option<f32>
    {
        Some(Into::<f32>::into(*self))
    }
    
    #[inline]
    fn to_f64(&self) -> Option<f64>
    {
        Some(Into::<f64>::into(*self))
    }
}