use std::ops::{Mul, MulAssign};

use crate::fp::{UInt, Fp, util};

use num_traits::NumCast;

impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize> Mul<Self> for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:
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

        if !e0.is_zero() || !Self::IS_INT_IMPLICIT
        {
            f0 = f0 + (self.int_bits() << FRAC_SIZE);
        }
        else
        {
            f0 = f0 << 1usize;
        }
        if !e1.is_zero() || !Self::IS_INT_IMPLICIT
        {
            f1 = f1 + (rhs.int_bits() << FRAC_SIZE);
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
                            f0 = util::rounding_div_2(f0);
                        }
                        else
                        {
                            f1 = util::rounding_div_2(f1);
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

        let mut o = U::zero();
        let mut f = loop
        {
            match f0.checked_mul(&f1)
            {
                Some(f) => break f,
                None => {
                    o = o + U::one();
                    if f0 > f1
                    {
                        f0 = util::rounding_div_2(f0);
                    }
                    else
                    {
                        f1 = util::rounding_div_2(f1);
                    }
                }
            }
        };

        let f_o = U::from(FRAC_SIZE - 1).unwrap();
        let mut e = match o.checked_sub(&f_o)
        {
            Some(o) => match e.checked_add(&o)
            {
                Some(e) => e,
                None => return if s {Self::neg_infinity()} else {Self::infinity()}
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
                None => return if s {Self::neg_infinity()} else {Self::infinity()}
            }
        };

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

        let s_bit = if s {U::one() << Self::SIGN_POS} else {U::zero()};

        if e.is_zero() && Self::IS_INT_IMPLICIT // subnormal
        {
            Fp::from_bits(s_bit + util::rounding_div_2(f))
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

            Fp::from_bits(s_bit + f + (e << Self::EXP_POS))
        }
    }
}
impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize> MulAssign for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:
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
        crate::tests::test_op2(Mul::mul, Mul::mul)
    }
}