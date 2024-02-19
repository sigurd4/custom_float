use std::ops::{Add, AddAssign};

use crate::{UInt, Fp, util};

impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize> Add<Self> for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:
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

            let e0 = self.exp_bits();
            let e1 = rhs.exp_bits();

            let mut f0 = self.frac_bits();
            let mut f1 = rhs.frac_bits();

            if !e0.is_zero() || !Self::IS_INT_IMPLICIT //normal
            {
                f0 = f0 + (self.int_bits() << FRAC_SIZE);
            }
            else
            {
                f0 = f0 << 1usize
            }
            if !e1.is_zero() || !Self::IS_INT_IMPLICIT //normal
            {
                f1 = f1 + (rhs.int_bits() << FRAC_SIZE);
            }
            else
            {
                f1 = f1 << 1usize
            }

            let mut e = if e0 >= e1
            {
                let shr = e0 - e1;
                f1 = util::rounding_div_pow(f1, U::from(2).unwrap(), shr);
                e0
            }
            else
            {
                let shr = e1 - e0;
                f0 = util::rounding_div_pow(f0, U::from(2).unwrap(), shr);
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
                            f0 = util::rounding_div_2(f0);
                            f1 = util::rounding_div_2(f1);
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

            while e > U::zero() && f <= U::one() << Self::MANTISSA_OP_SIZE - Self::BASE_SIZE
            {
                e = e - U::one();
                f = f << 1usize;
            }
            while f >= U::one() << Self::MANTISSA_OP_SIZE
            {
                e = e + U::one();
                f = util::rounding_div_2(f);
            }

            let n = if e.is_zero() && Self::IS_INT_IMPLICIT // subnormal
            {
                Fp::from_bits(util::rounding_div_2(f))
            }
            else
            {
                if Self::IS_INT_IMPLICIT
                {
                    f = f - (U::one() << FRAC_SIZE);
                }

                if e >= (U::one() << EXP_SIZE) - U::one()
                {
                    return if s {Self::neg_infinity()} else {Self::infinity()}
                }

                Fp::from_bits(f + (e << Self::EXP_POS))
            };
            if s {-n} else {n}
        }
    }
}
impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize> AddAssign for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:
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
        crate::tests::test_op2(Add::add, Add::add)
    }
}