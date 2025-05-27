use num_traits::ToBytes;

use crate::{util, Fp, UInt};

impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> ToBytes for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:,
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