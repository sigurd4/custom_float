use std::ops::{Div, DivAssign};

use crate::fp::{UInt, Fp, util};

use num_traits::{Inv, NumCast};

impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize> Div<Self> for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:
{
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output
    {
        if self.is_nan() || rhs.is_nan()
        {
            return self + rhs
        }
        if (self.is_zero() && rhs.is_zero()) || (self.is_infinite() && rhs.is_infinite())
        {
            return Self::qnan()
        }
    
        let s = !(self.sign_bit()^rhs.sign_bit()).is_zero();

        if self.is_infinite() || rhs.is_zero()
        {
            return if s {Self::neg_infinity()} else {Self::infinity()};
        }
        if rhs.is_infinite() || self.is_zero()
        {
            return if s {Self::neg_zero()} else {Self::zero()}
        }
        
        if rhs.is_one()
        {
            return self
        }
        if (-rhs).is_one()
        {
            return -self
        }

        let mut e0: U = self.exp_bits();
        let mut e1: U = rhs.exp_bits();

        let mut f0: U = self.frac_bits();
        let mut f1: U = rhs.frac_bits();

        if e0 == e1 && f0 == f1
        {
            return if s {-Self::one()} else {Self::one()}
        }

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

        let bias = Self::exp_bias();
        let e = match e0.checked_sub(&e1)
        {
            Some(e) => match e.checked_add(&bias)
            {
                Some(e) => e,
                None => return if s {Self::neg_infinity()} else {Self::infinity()}
            },
            None => {
                if e1 >= bias
                {
                    e1 = e1 - bias;
                }
                else
                {
                    e0 = e0 + (bias - e1);
                    e1 = U::zero();
                }
                match e0.checked_sub(&e1)
                {
                    Some(e) => e,
                    None => {
                        let mut o = e1 - e0;
                        while o > U::zero()
                        {
                            o = o - U::one();
                            if f1.leading_zeros() > 0
                            {
                                f1 = f1 << 1usize;
                            }
                            else
                            {
                                f0 = util::rounding_div_2(f0);
                            }
                        }
                        U::zero()
                    }
                }
            }
        };

        let mut o = U::from(f0.leading_zeros()).unwrap();
        f0 = f0 << f0.leading_zeros();
        let mut f = loop
        {
            if f1 == U::zero()
            {
                return if s {Self::neg_infinity()} else {Self::infinity()}
            }
            let f = f0 / f1;
            if !f0.is_zero() && f.leading_zeros() > 0
            {
                match f0.checked_shl(1)
                {
                    Some(f0_) => if (f0_ >> 1usize) + U::one() < f0
                    {
                        if f1.trailing_zeros() <= 0
                        {
                            break f
                        }
                        f1 = util::rounding_div_2(f1);
                    }
                    else
                    {
                        f0 = f0_
                    },
                    None => {
                        if f1.trailing_zeros() <= 0
                        {
                            break f
                        }
                        f1 = util::rounding_div_2(f1);
                    }
                }
                o = o + U::one();
            }
            else
            {
                break f
            }
        };

        let f_o = U::from(FRAC_SIZE).unwrap();
        let mut e = match f_o.checked_sub(&o)
        {
            Some(o) => match e.checked_add(&o)
            {
                Some(e) => e,
                None => return if s {Self::neg_infinity()} else {Self::infinity()}
            },
            None => match e.checked_add(&f_o)
            {
                Some(e) => match e.checked_sub(&o)
                {
                    Some(e) => e,
                    None => {
                        f = f >> <usize as NumCast>::from(o - e).unwrap();
                        U::zero()
                    }
                },
                None => return if s {Self::neg_infinity()} else {Self::infinity()}
            }
        };

        while e > U::zero() && f < U::one() << Self::MANTISSA_OP_SIZE - Self::BASE_SIZE
        {
            e = e - U::one();
            f = f << 1usize;
        }
        while e < U::one() << EXP_SIZE && f >= U::one() << Self::MANTISSA_OP_SIZE
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
impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize> DivAssign for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:
{
    #[inline]
    fn div_assign(&mut self, rhs: Self)
    {
        *self = *self / rhs
    }
}

impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize> Inv for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:
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
    use std::ops::Div;

    #[test]
    fn test_div()
    {
        crate::tests::test_op2(Div::div, Div::div)
    }
}