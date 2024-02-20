use num_traits::Pow;

use crate::{util, Fp, UInt};

macro_rules! impl_powu {
    ($($i:ty),*) => {
        $(
            impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Pow<$i> for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
            where
                [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
                [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:,
                [(); EXP_BASE - 2]:
            {
                type Output = Self;
            
                #[inline]
                fn pow(self, rhs: $i) -> Self::Output
                {
                    self.powu(rhs)
                }
            }
        )*
    };
}

impl_powu!(u8, u16, u32, usize, u64, u128/*, U256*/);

macro_rules! impl_powi {
    ($($i:ty),*) => {
        $(
            impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Pow<$i> for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
            where
                [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
                [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:,
                [(); EXP_BASE - 2]:
            {
                type Output = Self;
            
                #[inline]
                fn pow(self, rhs: $i) -> Self::Output
                {
                    self.powi(rhs)
                }
            }
        )*
    };
}

impl_powi!(i8, i16, i32, isize, i64, i128/*, I256*/);

impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Pow<Self> for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:,
    [(); EXP_BASE - 2]:
{
    type Output = Self;

    #[inline]
    fn pow(self, rhs: Self) -> Self::Output
    {
        self.powf(rhs)
    }
}