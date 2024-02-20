use crate::{util, Dec, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> PartialEq for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:
{
    fn eq(&self, other: &Self) -> bool
    {
        if self.is_nan() || other.is_nan()
        {
            return false
        }
        if self.is_zero() && other.is_zero()
        {
            return true;
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
        
        let ten = U::from(10).unwrap();
        while e0 > U::zero() && f0 <= U::one() << FRAC_SIZE - 1
        {
            e0 = e0 - U::one();
            f0 = f0*ten;
        }
        while f0 >= U::one() << FRAC_SIZE
        {
            e0 = e0 + U::one();
            f0 = f0/ten;
        }
        
        while e1 > U::zero() && f1 <= U::one() << FRAC_SIZE - 1
        {
            e1 = e1 - U::one();
            f1 = f1*ten;
        }
        while f1 >= U::one() << FRAC_SIZE
        {
            e1 = e1 + U::one();
            f1 = f1/ten;
        }

        e0 == e1 && f0 == f1
    }
}