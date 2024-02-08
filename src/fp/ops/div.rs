use std::ops::{Div, DivAssign, Mul};

use crate::fp::{UInt, Fp, bitsize_of};

use num_traits::{Float, Inv, NumCast, Zero};
use num_traits::One;

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> Div<Self> for Fp<U, EXP_SIZE, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 1]:
{
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output
    {
        if self.is_nan() || rhs.is_nan() || (self.is_zero() && rhs.is_zero()) || (self.is_infinite() && rhs.is_infinite())
        {
            return Self::nan()
        }
    
        let s = !(self.sign_bit()^rhs.sign_bit()).is_zero();

        if self.is_infinite() || rhs.is_zero()
        {
            return if s {Self::negative_infinity()} else {Self::infinity()};
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

        if !e0.is_zero() //normal
        {
            f0 = f0 + (U::one() << FRAC_SIZE);
        }
        else
        {
            f0 = f0 << 1usize;
        }
        if !e1.is_zero() //normal
        {
            f1 = f1 + (U::one() << FRAC_SIZE);
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
                None => return if s {Self::negative_infinity()} else {Self::infinity()}
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
                                f0 = f0 >> 1usize;
                            }
                            else
                            {
                                f1 = f1 >> 1usize;
                            }
                        }
                        U::zero()
                    }
                }
            }
        };

        f0 = f0 << bitsize_of::<U>() - FRAC_SIZE - 1;
        let mut o = U::from(bitsize_of::<U>() - FRAC_SIZE - 1).unwrap();
        let mut f = loop
        {
            if f1 == U::zero()
            {
                return if s {Self::negative_infinity()} else {Self::infinity()}
            }
            let f = f0 / f1;
            if !f0.is_zero() && f.leading_zeros() > 0
            {
                match f0.checked_shl(1)
                {
                    Some(f0_) => if (f0_ >> 1usize) + U::one() < f0
                    {
                        let f1_ = f1 >> 1usize;
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
                        let f1_ = f1 >> 1usize;
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

        let f_o = U::from(FRAC_SIZE).unwrap();
        let mut e = match f_o.checked_sub(&o)
        {
            Some(o) => match e.checked_add(&o)
            {
                Some(e) => e,
                None => return if s {Self::negative_infinity()} else {Self::infinity()}
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
                None => return if s {Self::negative_infinity()} else {Self::infinity()}
            }
        };

        while e < U::one() << EXP_SIZE && f >= U::one() << FRAC_SIZE + 1
        {
            e = e + U::one();
            f = f >> 1usize;
        }
        while e > U::zero() && f < U::one() << FRAC_SIZE
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
impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> DivAssign for Fp<U, EXP_SIZE, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 1]:
{
    fn div_assign(&mut self, rhs: Self)
    {
        *self = *self / rhs
    }
}

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> Inv for Fp<U, EXP_SIZE, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 1]:
{
    type Output = Self;

    fn inv(self) -> Self::Output
    {
        Self::one()/self
    }
}

#[cfg(test)]
mod test
{
    use std::ops::Div;

    use num_traits::{Float, One};

    use crate::fp::{ieee754::{FpDouble, FpSingle}, Fp};

    #[test]
    fn test_div()
    {
        //let n = FpSingle::from(3.333333)/FpSingle::from(1.0);
        //println!("{}", n);
        crate::tests::test_op2(Div::div, Div::div)
    }
}