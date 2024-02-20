use crate::{util, Dec, UInt};

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> PartialOrd for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:
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
        let mut e1 = self.exp_bits();
        let mut e2 = other.exp_bits();

        let mut f1 = self.frac_bits();
        let mut f2 = other.frac_bits();

        f1 = f1 + (self.int_bits() << FRAC_SIZE);
        f2 = f2 + (other.int_bits() << FRAC_SIZE);
        
        let ten = U::from(10).unwrap();
        while self.is_finite() && e1 > U::zero() && f1 <= U::one() << FRAC_SIZE - 3
        {
            e1 = e1 - U::one();
            f1 = f1*ten;
        }
        while self.is_finite() && e2 < ((U::one() << EXP_SIZE) - U::one() - U::one()) && f1 >= U::one() << FRAC_SIZE
        {
            e1 = e1 + U::one();
            f1 = util::rounding_div(f1, ten);
        }
        
        while other.is_finite() && e2 > U::zero() && f2 <= U::one() << FRAC_SIZE - 3
        {
            e2 = e2 - U::one();
            f2 = f2*ten;
        }
        while other.is_finite() && e2 < ((U::one() << EXP_SIZE) - U::one() - U::one()) && f2 >= U::one() << FRAC_SIZE
        {
            e2 = e2 + U::one();
            f2 = util::rounding_div(f2, ten);
        }
        
        if e1 != e2
        {
            return if s {e2.partial_cmp(&e1)} else {e1.partial_cmp(&e2)}
        }

        if s {f2.partial_cmp(&f1)} else {f1.partial_cmp(&f2)}
    }
}

#[cfg(test)]
mod test
{
    use num_traits::Float;

    #[test]
    fn test_rem()
    {
        crate::tests::test_op2_dec(Float::min, Float::min)
    }
}