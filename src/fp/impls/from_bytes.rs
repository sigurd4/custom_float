use num_traits::FromBytes;

use crate::{bitsize_of, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const INT_BIT: bool, const FRAC_SIZE: usize> FromBytes for Fp<U, EXP_SIZE, INT_BIT, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - INT_BIT as usize - FRAC_SIZE - 1]:,
    [(); bitsize_of::<U>() - EXP_SIZE - false as usize - FRAC_SIZE - 1]:,
    U: FromBytes
{
    type Bytes = U::Bytes;

    fn from_be_bytes(bytes: &Self::Bytes) -> Self
    {
        Self::from_bits(U::from_be_bytes(bytes))
    }

    fn from_le_bytes(bytes: &Self::Bytes) -> Self
    {
        Self::from_bits(U::from_le_bytes(bytes))
    }

    fn from_ne_bytes(bytes: &Self::Bytes) -> Self
    {
        Self::from_bits(U::from_ne_bytes(bytes))
    }
}