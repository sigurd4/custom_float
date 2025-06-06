use num_traits::FromBytes;

use crate::{Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> FromBytes for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE> + FromBytes
{
    type Bytes = U::Bytes;

    #[inline]
    fn from_be_bytes(bytes: &Self::Bytes) -> Self
    {
        Self::from_be_bytes(bytes)
    }

    #[inline]
    fn from_le_bytes(bytes: &Self::Bytes) -> Self
    {
        Self::from_le_bytes(bytes)
    }

    #[inline]
    fn from_ne_bytes(bytes: &Self::Bytes) -> Self
    {
        Self::from_ne_bytes(bytes)
    }
}