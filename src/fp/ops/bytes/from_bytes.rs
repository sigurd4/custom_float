use num_traits::FromBytes;

use crate::{util, Fp, UInt};

impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> FromBytes for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - 0 - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:,
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