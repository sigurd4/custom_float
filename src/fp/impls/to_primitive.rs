use num_traits::{NumCast, ToPrimitive};

use crate::{bitsize_of, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const INT_BIT: bool, const FRAC_SIZE: usize> ToPrimitive for Fp<U, EXP_SIZE, INT_BIT, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - INT_BIT as usize - FRAC_SIZE - 1]:,
    [(); bitsize_of::<U>() - EXP_SIZE - false as usize - FRAC_SIZE - 1]:
{
    fn to_i8(&self) -> Option<i8>
    {
        self.to_int()
    }
    fn to_u8(&self) -> Option<u8>
    {
        self.to_uint()
    }

    fn to_i16(&self) -> Option<i16>
    {
        self.to_int()
    }
    fn to_u16(&self) -> Option<u16>
    {
        self.to_uint()
    }

    fn to_i32(&self) -> Option<i32>
    {
        self.to_int()
    }
    fn to_u32(&self) -> Option<u32>
    {
        self.to_uint()
    }

    fn to_isize(&self) -> Option<isize>
    {
        self.to_int()
    }
    fn to_usize(&self) -> Option<usize>
    {
        self.to_uint()
    }

    fn to_i64(&self) -> Option<i64>
    {
        self.to_int()
    }
    fn to_u64(&self) -> Option<u64>
    {
        self.to_uint()
    }

    fn to_i128(&self) -> Option<i128>
    {
        self.to_int()
    }
    fn to_u128(&self) -> Option<u128>
    {
        self.to_uint()
    }

    fn to_f32(&self) -> Option<f32>
    {
        Some(Into::<f32>::into(*self))
    }

    fn to_f64(&self) -> Option<f64>
    {
        Some(Into::<f64>::into(*self))
    }
}