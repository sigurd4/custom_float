use core::num::FpCategory;

use crate::{Fp, FpRepr, UInt};

impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> PartialEq for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    fn eq(&self, other: &Self) -> bool
    {
        match (self.classify(), other.classify())
        {
            (FpCategory::Nan, _) | (_, FpCategory::Nan) => false,
            _ if self.0.eq(&other.0) => true,
            (FpCategory::Zero, FpCategory::Zero) => true,
            _ => {
                if Self::IS_INT_IMPLICIT
                {
                    return false
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
        
                Self::normalize_mantissa_down(&mut e0, &mut f0, Some(e1));
                Self::normalize_mantissa_down(&mut e1, &mut f1, Some(e0));
                Self::normalize_mantissa_up(&mut e0, &mut f0, Some(e1));
                Self::normalize_mantissa_up(&mut e1, &mut f1, Some(e0));
        
                e0 == e1 && f0 == f1
            }
        }
    }
}