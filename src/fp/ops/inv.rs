use num_traits::Inv;

use crate::{util, Fp, UInt};

impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Inv for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:
{
    type Output = Self;

    #[inline]
    fn inv(self) -> Self::Output
    {
        self.recip()
    }
}

#[cfg(test)]
mod test
{
    use num_traits::Inv;

    #[test]
    fn test_inv()
    {
        crate::tests::test_op1("inv", Inv::inv, Inv::inv, Some(0.00001), Some(0.1..10.0))
    }
}