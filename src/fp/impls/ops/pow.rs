use num_traits::Pow;

use crate::{bitsize_of, Fp, UInt};

macro_rules! impl_powu {
    ($($i:ty),*) => {
        $(
            impl<U: UInt, const EXP_SIZE: usize, const INT_BIT: bool, const FRAC_SIZE: usize> Pow<$i> for Fp<U, EXP_SIZE, INT_BIT, FRAC_SIZE>
            where
                [(); bitsize_of::<U>() - EXP_SIZE - INT_BIT as usize - FRAC_SIZE - 1]:,
                [(); bitsize_of::<U>() - EXP_SIZE - false as usize - FRAC_SIZE - 1]:
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
            impl<U: UInt, const EXP_SIZE: usize, const INT_BIT: bool, const FRAC_SIZE: usize> Pow<$i> for Fp<U, EXP_SIZE, INT_BIT, FRAC_SIZE>
            where
                [(); bitsize_of::<U>() - EXP_SIZE - INT_BIT as usize - FRAC_SIZE - 1]:,
                [(); bitsize_of::<U>() - EXP_SIZE - false as usize - FRAC_SIZE - 1]:
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

impl<U: UInt, const EXP_SIZE: usize, const INT_BIT: bool, const FRAC_SIZE: usize> Pow<Self> for Fp<U, EXP_SIZE, INT_BIT, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - INT_BIT as usize - FRAC_SIZE - 1]:,
    [(); bitsize_of::<U>() - EXP_SIZE - false as usize - FRAC_SIZE - 1]:
{
    type Output = Self;

    #[inline]
    fn pow(self, rhs: Self) -> Self::Output
    {
        self.powf(rhs)
    }
}