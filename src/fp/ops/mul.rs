use std::ops::{Mul, MulAssign};

use crate::fp::{UInt, Fp, bitsize_of};

use num_traits::{Float, NumCast, One, Zero};

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> Mul<Self> for Fp<U, EXP_SIZE, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 1]:
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output
    {
        if self.is_nan() || rhs.is_nan() || (self.is_infinite() && rhs.is_zero()) || (self.is_zero() && rhs.is_infinite())
        {
            return Self::nan()
        }
        
        if self.is_one()
        {
            return rhs
        }
        if (-self).is_one()
        {
            return -rhs
        }
        if rhs.is_one()
        {
            return self
        }
        if (-rhs).is_one()
        {
            return -self
        }
    
        let s = !(self.sign_bit()^rhs.sign_bit()).is_zero();

        if self.is_zero() || rhs.is_zero()
        {
            return if s {-Self::zero()} else {Self::zero()}
        }
        if self.is_infinite() || rhs.is_infinite()
        {
            return if s {Self::negative_infinity()} else {Self::infinity()}
        }
    
        let mut e0: U = self.exp_bits();
        let mut e1: U = rhs.exp_bits();

        let mut f0: U = self.frac_bits();
        let mut f1: U = rhs.frac_bits();

        if !e0.is_zero()
        {
            f0 = f0 + (U::one() << FRAC_SIZE);
        }
        else
        {
            f0 = f0 << 1usize;
        }
        if !e1.is_zero()
        {
            f1 = f1 + (U::one() << FRAC_SIZE);
        }
        else
        {
            f1 = f1 << 1usize;
        }

        let bias = Self::exp_bias() + U::one();
        let e = match e0.checked_add(&e1)
        {
            Some(e) => match e.checked_sub(&bias)
            {
                Some(e) => e,
                None => {
                    let mut o = bias - e;
                    while o > U::zero()
                    {
                        o = o - U::one();
                        if f0 > f1
                        {
                            f0 = f0 >> 1usize;
                        }
                        else
                        {
                            f1 = f1 >> 1usize;
                        }
                    }
                    U::zero()
                }
            },
            None => {
                if e0 > e1
                {
                    e0 = e0 - bias;
                }
                else
                {
                    e1 = e1 - bias;
                }
                match e0.checked_add(&e1)
                {
                    Some(e) => e,
                    None => return if s {Self::negative_infinity()} else {Self::infinity()}
                }
            }
        };

        let mut o = U::zero();
        let mut f = loop
        {
            let (f, c) = f0.overflowing_mul(&f1);
            if c
            {
                o = o + U::one();
                if f0 > f1
                {
                    f0 = f0 >> 1usize;
                }
                else
                {
                    f1 = f1 >> 1usize;
                }
            }
            else
            {
                break f
            }
        };

        let f_o = U::from(FRAC_SIZE - 1).unwrap();
        let mut e = match o.checked_sub(&f_o)
        {
            Some(o) => match e.checked_add(&o)
            {
                Some(e) => e,
                None => return if s {Self::negative_infinity()} else {Self::infinity()}
            },
            None => match e.checked_add(&o)
            {
                Some(e) => match e.checked_sub(&f_o)
                {
                    Some(e) => e,
                    None => {
                        f = f >> <usize as NumCast>::from(f_o - e).unwrap();
                        U::zero()
                    }
                },
                None => return if s {Self::negative_infinity()} else {Self::infinity()}
            }
        };

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

        let s_bit = if s {U::one() << Self::SIGN_POS} else {U::zero()};

        if e.is_zero() // subnormal
        {
            return Fp::from_bits(s_bit + (f >> 1usize))
        }
        else
        {
            if e >= (U::one() << EXP_SIZE) - U::one()
            {
                return if s {Self::negative_infinity()} else {Self::infinity()}
            }

            f = f - (U::one() << FRAC_SIZE);

            return Fp::from_bits(s_bit + f + (e << Self::EXP_POS))
        }
    }
}
impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> MulAssign for Fp<U, EXP_SIZE, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 1]:
{
    fn mul_assign(&mut self, rhs: Self)
    {
        *self = *self * rhs
    }
}

#[cfg(test)]
mod test
{
    use std::ops::Mul;

    use crate::fp::ieee754::FpSingle;

    #[test]
    fn test_mul()
    {
        crate::tests::test_op2(Mul::mul, Mul::mul)
    }
}