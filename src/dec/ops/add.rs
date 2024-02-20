use std::ops::{Add, AddAssign};

use num_traits::NumCast;

use crate::{util, Dec, UInt};

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> Add<Self> for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output
    {
        if self.is_nan() || rhs.is_nan()
        {
            if rhs.is_snan() || !self.is_nan()
            {
                return rhs
            }
            self
        }
        else if self.is_infinite()
        {
            if rhs.is_infinite()
            {
                if self.sign_bit() == rhs.sign_bit()
                {
                    self
                }
                else
                {
                    Self::qnan()
                }
            }
            else
            {
                self
            }
        }
        else if rhs.is_infinite()
        {
            rhs
        }
        else if self.is_zero()
        {
            rhs
        }
        else if rhs.is_zero()
        {
            self
        }
        else
        {
            let s0 = !self.sign_bit().is_zero();
            let s1 = !rhs.sign_bit().is_zero();

            let mut e0 = self.exp_bits();
            let mut e1 = rhs.exp_bits();

            let mut f0 = self.frac_bits();
            let mut f1 = rhs.frac_bits();

            f0 = f0 + (self.int_bits() << FRAC_SIZE);
            f1 = f1 + (rhs.int_bits() << FRAC_SIZE);

            let ten = U::from(10).unwrap();

            let mut e = if e0 >= e1
            {
                while e1 < e0
                {
                    f1 = util::rounding_div(f1, ten);
                    e1 = e1 + U::one();
                }
                e0
            }
            else
            {
                while e0 < e1
                {
                    f0 = util::rounding_div(f0, ten);
                    e0 = e0 + U::one();
                }
                e1
            };

            let s = match (s0, s1)
            {
                (false, false) => false,
                (false, true) => f0 < f1,
                (true, false) => f0 > f1,
                (true, true) => true,
            };

            let mut f = if s0 == s1
            {
                loop
                {
                    match f0.checked_add(&f1)
                    {
                        Some(f) => break f,
                        None => {
                            e = e + U::one();
                            f0 = util::rounding_div(f0, ten);
                            f1 = util::rounding_div(f1, ten);
                        }
                    }
                }
            }
            else
            {
                if f0 >= f1 {f0 - f1} else {f1 - f0}
            };

            if f.is_zero()
            {
                return if s {-Self::zero()} else {Self::zero()}
            }

            while e > U::zero() && f <= U::one() << FRAC_SIZE - 4
            {
                e = e - U::one();
                f = f*ten;
            }
            while f >= U::one() << FRAC_SIZE
            {
                e = e + U::one();
                f = util::rounding_div(f, ten)
            }

            if e >= (U::one() << EXP_SIZE) - U::one()
            {
                return if s {Self::neg_infinity()} else {Self::infinity()}
            }

            let n = Self::from_bits(f + (e << Self::EXP_POS));

            if s {-n} else {n}
        }
    }
}
impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> AddAssign for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:
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

    #[test]
    fn test_add()
    {
        crate::tests::test_op2_dec(Add::add, Add::add)
    }
}