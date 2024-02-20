use num_traits::NumCast;

use crate::{util, Dec, UInt};

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> NumCast for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:
{
    #[inline]
    fn from<T: num_traits::ToPrimitive>(n: T) -> Option<Self>
    {
        n.to_f64().map(|f| f.into())
    }
}