use std::ops::{Mul, MulAssign};

use crate::{util, Fp, Dec, UInt};

use num_traits::NumCast;

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> Mul<Self> for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output
    {
        if self.is_nan() || rhs.is_nan()
        {
            return self + rhs
        }
        if (self.is_infinite() && rhs.is_zero()) || (self.is_zero() && rhs.is_infinite())
        {
            return Self::qnan()
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
            return if s {Self::neg_infinity()} else {Self::infinity()}
        }
    
        let mut e0: U = self.exp_bits();
        let mut e1: U = rhs.exp_bits();

        let mut f0: U = self.frac_bits();
        let mut f1: U = rhs.frac_bits();

        f0 = f0 + (self.int_bits() << Self::INT_POS);
        f1 = f1 + (rhs.int_bits() << Self::INT_POS);

        let ten = U::from(10).unwrap();
        let bias = Self::exp_bias() + U::one();
        let mut e = match e0.checked_add(&e1)
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
                            f0 = util::rounding_div(f0, ten);
                        }
                        else
                        {
                            f1 = util::rounding_div(f1, ten);
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
                    None => return if s {Self::neg_infinity()} else {Self::infinity()}
                }
            }
        };

        let mut o = U::one();
        let mut f = loop
        {
            match f0.checked_mul(&f1)
            {
                Some(f) => break f,
                None => {
                    o = o + U::one();
                    if f0 > f1
                    {
                        f0 = util::rounding_div(f0, ten);
                    }
                    else
                    {
                        f1 = util::rounding_div(f1, ten);
                    }
                }
            }
        };

        for _ in 0..FRAC_SIZE
        {
            while (o > U::zero() || e > U::zero()) && f.leading_zeros() > 4
            {
                if o > U::zero()
                {
                    o = o - U::one();
                }
                else
                {
                    e = e - U::one();
                }
                f = f*ten;
            }
            f = f >> 1usize;
        }
        
        let mut e = match e.checked_add(&o)
        {
            Some(e) => e,
            None => return if s {Self::neg_infinity()} else {Self::infinity()}
        };

        while e > U::zero() && f <= U::one() << FRAC_SIZE - 4
        {
            e = e - U::one();
            f = f*ten;
        }
        while f >= U::one() << FRAC_SIZE
        {
            e = e + U::one();
            f = util::rounding_div(f, ten);
        }

        let s_bit = if s {U::one() << Self::SIGN_POS} else {U::zero()};

        if e >= (U::one() << EXP_SIZE) - U::one()
        {
            return if s {Self::neg_infinity()} else {Self::infinity()}
        }

        Self::from_bits(s_bit + f + (e << Self::EXP_POS))
    }
}
impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> MulAssign for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:
{
    #[inline]
    fn mul_assign(&mut self, rhs: Self)
    {
        *self = *self * rhs
    }
}

#[cfg(test)]
mod test
{
    use std::ops::Mul;

    #[test]
    fn test_mul()
    {
        crate::tests::test_op2_dec(Mul::mul, Mul::mul)
    }
}