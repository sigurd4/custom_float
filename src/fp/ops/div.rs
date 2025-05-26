use core::ops::{Div, DivAssign};

use crate::fp::{UInt, Fp, util};

impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Div<Self> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - 0 - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:
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
    
        let s = self.is_sign_negative()^rhs.is_sign_negative();

        if self.is_infinite() || rhs.is_zero()
        {
            return if s {Self::neg_infinity()} else {Self::infinity()};
        }
        if rhs.is_infinite() || self.is_zero()
        {
            return if s {Self::neg_zero()} else {Self::zero()}
        }
        
        if rhs.abs().is_one()
        {
            return self.with_sign(s)
        }

        let mut e0: U = self.exp_bits();
        let mut e1: U = rhs.exp_bits();

        let mut f0: U = self.mantissa_bits();
        let mut f1: U = rhs.mantissa_bits();

        if e0 == e1 && f0 == f1
        {
            return if s {-Self::one()} else {Self::one()}
        }

        let base = U::from(EXP_BASE).unwrap();
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
                            if (f0 % base).is_zero()
                            {
                                f0 = f0/base;
                            }
                            else if f1.leading_zeros() as usize > Self::BASE_PADDING
                            {
                                f1 = f1*base;
                            }
                            else
                            {
                                f0 = util::rounding_div(f0, base);
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
            if f1.is_zero()
            {
                return if s {Self::neg_infinity()} else {Self::infinity()}
            }
            let f = util::rounding_div(f0, f1);
            if !f0.is_zero() && f.leading_zeros() as usize >= Self::BASE_PADDING
            {
                if f0.leading_zeros() as usize > Self::BASE_PADDING
                {
                    f0 = f0*base
                }
                else
                {
                    if f1 % base != U::zero()
                    {
                        break f
                    }
                    f1 = f1/base;
                }
                o = o + U::one();
            }
            else
            {
                break f
            }
        };

        if EXP_BASE == 2
        {
            let frac_size = U::from(FRAC_SIZE).unwrap();
            if o >= frac_size
            {
                o = o - frac_size
            }
            else{
                match e.checked_add(&frac_size)
                {
                    Some(e_) => e = e_,
                    None => return if s {Self::neg_infinity()} else {Self::infinity()}
                }
            }
        }
        else
        {
            for _ in 0..FRAC_SIZE
            {
                f = loop
                {
                    if f.leading_zeros() > 0
                    {
                        break f << 1usize
                    }
                
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
                    f = util::rounding_div(f, base);
                }
            }
        }

        while e < U::one() << EXP_SIZE && f >= U::one() << Self::MANTISSA_OP_SIZE
        {
            e = e + U::one();
            f = util::rounding_div(f, base);
        }

        let mut e = match e.checked_sub(&o)
        {
            Some(e) => e,
            None => return if s {Self::neg_zero()} else {Self::zero()}
        };

        Self::carry_exp_mantissa(&mut e, &mut f);
        Self::from_sign_exp_mantissa(s, e, f)
    }
}
impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> DivAssign for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - 0 - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:
{
    #[inline]
    fn div_assign(&mut self, rhs: Self)
    {
        *self = *self / rhs
    }
}

#[cfg(test)]
mod test
{
    use std::ops::Div;

    use crate::tests::F;

    #[test]
    fn test_div_once()
    {
        let a = F::from(-2.2e0f32);
        let b = F::from(3.333333e0f32);
        let c = a / b;
        println!("{a} / {b} = {c}");
    }

    #[test]
    fn test_div()
    {
        crate::tests::test_op2("div", Div::div, Div::div, None)
    }
}