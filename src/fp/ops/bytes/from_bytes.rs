use num_traits::FromBytes;

use crate::{util, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize> FromBytes for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:,
    U: FromBytes
{
    type Bytes = U::Bytes;

    #[inline]
    fn from_be_bytes(bytes: &Self::Bytes) -> Self
    {
        Self::from_bits(U::from_be_bytes(bytes))
    }

    #[inline]
    fn from_le_bytes(bytes: &Self::Bytes) -> Self
    {
        Self::from_bits(U::from_le_bytes(bytes))
    }

    #[inline]
    fn from_ne_bytes(bytes: &Self::Bytes) -> Self
    {
        Self::from_bits(U::from_ne_bytes(bytes))
    }
}