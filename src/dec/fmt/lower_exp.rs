use std::fmt::LowerExp;

use crate::{util, Dec, UInt};

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> LowerExp for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result
    {
        let x = Into::<f64>::into(*self);
        if self.is_finite() != x.is_finite()
        {
            write!(f, "~")?
        }
        <f64 as LowerExp>::fmt(&x, f)
    }
}