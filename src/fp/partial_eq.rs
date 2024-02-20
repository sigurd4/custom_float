use crate::{util, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> PartialEq for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:,
    [(); EXP_BASE - 2]:
{
    fn eq(&self, other: &Self) -> bool
    {
        if self.is_nan() || other.is_nan()
        {
            return false
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
        
        while e0 > e1 && f0 <= U::one() << Self::MANTISSA_OP_SIZE - Self::BASE_SIZE
        {
            e0 = e0 - U::one();
            f0 = f0 << 1usize;
        }
        while e0 < e1 && f0 >= U::one() << Self::MANTISSA_OP_SIZE
        {
            e0 = e0 + U::one();
            f0 = f0 >> 1usize;
        }
        
        while e1 > e0 && f1 <= U::one() << Self::MANTISSA_OP_SIZE - Self::BASE_SIZE
        {
            e1 = e1 - U::one();
            f1 = f1 << 1usize;
        }
        while e1 < e0 && f1 >= U::one() << Self::MANTISSA_OP_SIZE
        {
            e1 = e1 + U::one();
            f1 = f1 >> 1usize;
        }

        e0 == e1 && f0 == f1
    }
}