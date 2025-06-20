use core::ops::{Add, AddAssign};

use crate::{Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Add<Self> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output
    {
        self.add_with_sign(false, rhs, false)
    }
}
impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> AddAssign for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
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
    use core::ops::Add;
    use test::Bencher;

    use crate::{ieee754::FpDouble, tests::{bench_op2, test_op2}};

    #[test]
    fn test_add_once()
    {
        type F = FpDouble;

        let a = F::from(1f32);
        let b = F::from(-16f32);
        let c = a + b;
        println!("{a} + {b} = {c}");
    }
    #[test]
    fn test_add()
    {
        test_op2!("add", Add::add, Some(0.001))
    }
    #[bench]
    fn bench_add(bencher: &mut Bencher)
    {
        test_add();
        bench_op2!(bencher, Add::add)
    }
}