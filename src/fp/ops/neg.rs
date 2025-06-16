use core::ops::Neg;

use crate::{Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Neg for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output
    {
        if !SIGN_BIT
        {
            if self.is_zero()
            {
                return self
            }
            return Self::qnan()
        }
        Self::from_bits(self.to_bits().bitxor(Self::shift_sign(U::one())))
    }
}

#[cfg(test)]
mod test
{
    use core::ops::Neg;
    use test::Bencher;

    use crate::tests::{bench_op1, test_op1};

    #[test]
    fn test_neg()
    {
        test_op1!("neg", Neg::neg, None, Some(-10.0..10.0))
    }
    #[bench]
    fn bench_neg(bencher: &mut Bencher)
    {
        test_neg();
        bench_op1!(bencher, Neg::neg)
    }
}