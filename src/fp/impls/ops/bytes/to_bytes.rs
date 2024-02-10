use num_traits::ToBytes;

use crate::{bitsize_of, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const INT_BIT: bool, const FRAC_SIZE: usize> ToBytes for Fp<U, EXP_SIZE, INT_BIT, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - INT_BIT as usize - FRAC_SIZE - 1]:,
    [(); bitsize_of::<U>() - EXP_SIZE - false as usize - FRAC_SIZE - 1]:,
    U: ToBytes
{
    type Bytes = U::Bytes;

    #[inline]
    fn to_be_bytes(&self) -> Self::Bytes
    {
        self.to_bits().to_be_bytes()
    }

    #[inline]
    fn to_le_bytes(&self) -> Self::Bytes
    {
        self.to_bits().to_le_bytes()
    }

    #[inline]
    fn to_ne_bytes(&self) -> Self::Bytes
    {
        self.to_bits().to_ne_bytes()
    }
}