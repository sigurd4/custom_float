use num_traits::NumCast;

use crate::{bitsize_of, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> NumCast for Fp<U, EXP_SIZE, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 1]:
{
    fn from<T: num_traits::ToPrimitive>(n: T) -> Option<Self>
    {
        n.to_f64()
            .map(|f| f.into())
    }
}