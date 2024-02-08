use num_traits::Float;

use crate::{bitsize_of, Fp, UInt};


impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> PartialOrd for Fp<U, EXP_SIZE, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 1]:
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering>
    {
        if self.is_nan() || other.is_nan()
        {
            return None
        }
        let s1 = self.sign_bit();
        let s2 = other.sign_bit();
        if s1 != s2
        {
            return s2.partial_cmp(&s1)
        }
        let s = !s1.is_zero();
        let e1 = self.exp_bits();
        let e2 = other.exp_bits();
        if e1 != e2
        {
            return if s {e2.partial_cmp(&e1)} else {e1.partial_cmp(&e2)}
        }
        let f1 = self.frac_bits();
        let f2 = other.frac_bits();
        return if s {f2.partial_cmp(&f1)} else {f1.partial_cmp(&f2)}
    }
}