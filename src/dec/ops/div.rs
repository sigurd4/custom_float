use std::ops::{Div, DivAssign};

use crate::{util, Fp, UInt, Dec};

use num_traits::{Euclid, Inv, NumCast};

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> Div<Self> for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:
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

        f0 = f0 + (self.int_bits() << FRAC_SIZE);
        f1 = f1 + (rhs.int_bits() << FRAC_SIZE);

        let ten = U::from(10).unwrap();
        let bias = Self::exp_bias();
        let mut e = match e0.checked_sub(&e1)
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
                }
            }
        };

        let mut o = U::zero();
        let mut f = loop
        {
            if f1 == U::zero()
            {
                return if s {Self::neg_infinity()} else {Self::infinity()}
            }
            let f = f0 / f1;
            if !f0.is_zero() && f.leading_zeros() > 0
            {
                match f0.checked_mul(&ten)
                {
                    Some(f0_) => if util::rounding_div(f0_, ten) + U::one() < f0
                    {
                        let f1_ = util::rounding_div(f1, ten);
                        if f1.trailing_zeros() <= 0
                        {
                            break f
                        }
                        else
                        {
                            f1 = f1_;
                        }
                    }
                    else
                    {
                        f0 = f0_
                    },
                    None => {
                        let f1_ = util::rounding_div(f1, ten);
                        if f1.trailing_zeros() <= 0
                        {
                            break f
                        }
                        else
                        {
                            f1 = f1_;
                        }
                    }
                }
                o = o + U::one();
            }
            else
            {
                break f
            }
        };

        for _ in 0..FRAC_SIZE
        {
            f = loop
            {
                match f.checked_shl(1)
                {
                    Some(f_) => if (f_ >> 1usize) + U::one() < f
                    {
                        if f.is_zero()
                        {
                            return Self::zero()
                        }
                        if o > U::zero()
                        {
                            o = o - U::one();
                        }
                        else if e < U::max_value()
                        {
                            e = e + U::one();
                        }
                        else
                        {
                            return if s {Self::neg_infinity()} else {Self::infinity()}
                        }
                        f = util::rounding_div(f, ten);
                    }
                    else
                    {
                        break f_;
                    },
                    None => {
                        if f.is_zero()
                        {
                            return Self::zero()
                        }
                        if o > U::zero()
                        {
                            o = o - U::one();
                        }
                        else if e < U::max_value()
                        {
                            e = e + U::one();
                        }
                        else
                        {
                            return if s {Self::neg_infinity()} else {Self::infinity()}
                        }
                        f = util::rounding_div(f, ten);
                    }
                }
            }
        }

        while e < U::one() << EXP_SIZE && f >= U::one() << FRAC_SIZE
        {
            e = e + U::one();
            f = util::rounding_div(f, ten);
        }

        let mut e = match e.checked_sub(&o)
        {
            Some(e) => e,
            None => return if s {Self::neg_zero()} else {Self::zero()}
        };

        while e > U::zero() && f < U::one() << FRAC_SIZE - 3
        {
            e = e - U::one();
            f = f*ten;
        }
        while e < U::one() << EXP_SIZE && f >= U::one() << FRAC_SIZE
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
impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> DivAssign for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:
{
    #[inline]
    fn div_assign(&mut self, rhs: Self)
    {
        *self = *self / rhs
    }
}

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> Inv for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:
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

    use num_traits::ToPrimitive;

    use crate::ieee754::DecDouble;

    #[test]
    fn test_div()
    {
        crate::tests::test_op2_dec(Div::div, Div::div)
    }
}