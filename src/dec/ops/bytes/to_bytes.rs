use num_traits::ToBytes;

use crate::{util, Dec, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> ToBytes for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:,
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