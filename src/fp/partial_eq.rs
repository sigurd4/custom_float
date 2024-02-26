use crate::{util, Fp, UInt};

impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> PartialEq for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - 0 - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:
{
    fn eq(&self, other: &Self) -> bool
    {
        if self.is_nan() || other.is_nan()
        {
            return false
        }
        if self.is_zero() && other.is_zero()
        {
            return true
        }
        if Self::IS_INT_IMPLICIT
        {
            return self.0.eq(&other.0)
        }

        let s0 = self.sign_bit();
        let s1 = other.sign_bit();

        if s0 != s1
        {
            return false;
        }
    
        let mut e0: U = self.exp_bits();
        let mut e1: U = other.exp_bits();

        let mut f0: U = self.frac_bits();
        let mut f1: U = other.frac_bits();

        f0 = f0 + (self.int_bits() << FRAC_SIZE);
        f1 = f1 + (other.int_bits() << FRAC_SIZE);

        let base = U::from(EXP_BASE).unwrap();
        
        while e0 > e1 && f0 <= U::one() << Self::MANTISSA_OP_SIZE - Self::BASE_PADDING
        {
            e0 = e0 - U::one();
            f0 = f0*base;
        }
        while e1 > e0 && f1 <= U::one() << Self::MANTISSA_OP_SIZE - Self::BASE_PADDING
        {
            e1 = e1 - U::one();
            f1 = f1*base;
        }

        while e0 < e1 && f0 >= U::one() << Self::MANTISSA_OP_SIZE
        {
            e0 = e0 + U::one();
            f0 = util::rounding_div(f0, base);
        }
        while e1 < e0 && f1 >= U::one() << Self::MANTISSA_OP_SIZE
        {
            e1 = e1 + U::one();
            f1 = util::rounding_div(f1, base);
        }

        e0 == e1 && f0 == f1
    }
}