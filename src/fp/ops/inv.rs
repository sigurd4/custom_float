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
    use test::Bencher;

    use crate::tests::{self, F};

    #[test]
    fn test_inv()
    {
        crate::tests::test_op1("inv", Inv::inv, Inv::inv, None, Some(0.1..10.0))
    }
    #[bench]
    fn bench_inv(bencher: &mut Bencher)
    {
        test_inv();
        tests::bench_op1::<F, _>(bencher, Inv::inv)
    }
}