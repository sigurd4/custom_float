use core::ops::{Mul, MulAssign};

use crate::fp::{UInt, Fp, util};

impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Mul<Self> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - 0 - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:
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

        let s = self.is_sign_negative()^rhs.is_sign_negative();
        
        if self.abs().is_one()
        {
            return rhs.with_sign(s)
        }
        if rhs.abs().is_one()
        {
            return self.with_sign(s)
        }

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

        let mut f0: U = self.mantissa_bits();
        let mut f1: U = rhs.mantissa_bits();

        let base = U::from(EXP_BASE).unwrap();
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
                        if (f0 % base).is_zero()
                        {
                            f0 = f0/base
                        }
                        else if (f1 % base).is_zero()
                        {
                            f1 = f1/base
                        }
                        else if f0 > f1
                        {
                            f0 = util::rounding_div(f0, base);
                        }
                        else
                        {
                            f1 = util::rounding_div(f1, base);
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
                    if (f0 % base).is_zero()
                    {
                        f0 = f0/base
                    }
                    else if (f1 % base).is_zero()
                    {
                        f1 = f1/base
                    }
                    else if f0 > f1
                    {
                        f0 = util::rounding_div(f0, base);
                    }
                    else
                    {
                        f1 = util::rounding_div(f1, base);
                    }
                }
            }
        };

        for _ in 0..FRAC_SIZE
        {
            while (o > U::zero() || e > U::zero()) && f.leading_zeros() as usize > Self::BASE_PADDING
            {
                if o > U::zero()
                {
                    o = o - U::one();
                }
                else
                {
                    e = e - U::one();
                }
                f = f*base;
            }
            f = f >> 1usize;
        }
        
        let mut e = match e.checked_add(&o)
        {
            Some(e) => e,
            None => return if s {Self::neg_infinity()} else {Self::infinity()}
        };

        Self::carry_exp_mantissa(&mut e, &mut f);
        Self::from_sign_exp_mantissa(s, e, f)
    }
}
impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> MulAssign for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - 0 - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:
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
        crate::tests::test_op2("mul", Mul::mul, Mul::mul, None)
    }
}