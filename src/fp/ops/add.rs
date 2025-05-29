use core::ops::{Add, AddAssign};

use crate::{UInt, Fp, util};

impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Add<Self> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:
{
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output
    {
        self.add_with_sign(rhs, false)
    }
}
impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> AddAssign for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:
{
    #[inline]
    fn add_assign(&mut self, rhs: Self)
    {
        *self = *self + rhs
    }
}

#[cfg(test)]
mod test
{
    use std::ops::Add;

    use test::Bencher;

    use crate::tests::{self, F};

    #[test]
    fn test_add_once()
    {
        let a = F::from(-2.2e0f32);
        let b = F::from(3.333333e0f32);
        let c = a + b;
        println!("{a} + {b} = {c}");
    }
    #[test]
    fn test_add()
    {
        crate::tests::test_op2("add", Add::add, Add::add, None)
    }
    #[bench]
    fn bench_add(bencher: &mut Bencher)
    {
        test_add();
        tests::bench_op2::<F, _>(bencher, Add::add)
    }
}