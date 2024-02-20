use num_traits::FromBytes;

use crate::{util, Dec, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> FromBytes for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:,
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