use num_traits::Euclid;

use crate::{Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Euclid for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[inline]
    fn div_euclid(&self, rhs: &Self) -> Self
    {
        (*self).div_euclid(*rhs)
    }

    #[inline]
    fn rem_euclid(&self, rhs: &Self) -> Self
    {
        (*self).rem_euclid(*rhs)
    }
}