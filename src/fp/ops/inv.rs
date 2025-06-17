use num_traits::Inv;

use crate::{Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Inv for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
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

    use crate::{khronos::KhronosFp10, tests::{bench_op1, test_op1}};

    #[test]
    fn test_inv_once()
    {
        type F = KhronosFp10;

        let x = -F::one();
        let y = x.inv();

        println!("1 / {x} = {y}")
    }

    #[test]
    fn test_inv()
    {
        test_op1!("inv", Inv::inv, None, Some(0.1..10.0))
    }
    #[bench]
    fn bench_inv(bencher: &mut Bencher)
    {
        test_inv();
        bench_op1!(bencher, Inv::inv)
    }
}