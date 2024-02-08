use std::ops::{Add, AddAssign};

use num_traits::{Float, NumCast, Zero};

use crate::fp::ieee754::FpSingle;

use super::super::{UInt, Fp, bitsize_of};

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> Add<Self> for Fp<U, EXP_SIZE, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 1]:
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output
    {
        if self.is_nan() || rhs.is_nan()
        {
            Self::nan()
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
                    Self::nan()
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

            if !e0.is_zero() //normal
            {
                f0 = f0 + (U::one() << FRAC_SIZE);
            }
            else
            {
                f0 = f0 << 1usize
            }
            if !e1.is_zero() //normal
            {
                f1 = f1 + (U::one() << FRAC_SIZE);
            }
            else
            {
                f1 = f1 << 1usize
            }

            let mut e = if e0 >= e1
            {
                let shr = <u32 as NumCast>::from(e0 - e1).unwrap();
                f1 = match f1.checked_shr(shr)
                {
                    Some(f1) => f1,
                    None => U::zero()
                };
                e0
            }
            else
            {
                let shr = <u32 as NumCast>::from(e1 - e0).unwrap();
                f0 = match f0.checked_shr(shr)
                {
                    Some(f0) => f0,
                    None => U::zero()
                };
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
                            f0 = f0 >> 1usize;
                            f1 = f1 >> 1usize;
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

            while f >= U::one() << FRAC_SIZE + 1
            {
                e = e + U::one();
                f = f >> 1usize;
            }
            while e > U::zero() && f <= U::one() << FRAC_SIZE
            {
                e = e - U::one();
                f = f << 1usize;
            }

            let n = if e.is_zero() // subnormal
            {
                Fp::from_bits(f >> 1usize)
            }
            else
            {
                if e >= (U::one() << EXP_SIZE) - U::one()
                {
                    return if s {Self::negative_infinity()} else {Self::infinity()}
                }

                f = f - (U::one() << FRAC_SIZE);

                Fp::from_bits(f + (e << FRAC_SIZE))
            };
            if s {-n} else {n}
        }
    }
}
impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> AddAssign for Fp<U, EXP_SIZE, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 1]:
{
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