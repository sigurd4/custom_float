use core::ops::{Sub, SubAssign};

use crate::{Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Sub<Self> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output
    {
        self.add_with_sign(false, rhs, true)
    }
}
impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> SubAssign for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[inline]
    fn sub_assign(&mut self, rhs: Self)
    {
        *self = *self - rhs
    }
}

#[cfg(test)]
mod test
{
    use core::ops::Sub;
    use test::Bencher;

    use crate::tests::{bench_op2, test_op2};

    #[test]
    fn test_sub()
    {
        test_op2!("sub", Sub::sub, None)
    }
    #[bench]
    fn bench_sub(bencher: &mut Bencher)
    {
        test_sub();
        bench_op2!(bencher, Sub::sub)
    }
}