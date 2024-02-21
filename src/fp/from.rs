use crate::{ieee754::{FpDouble, FpSingle}, util, Fp, UInt};

macro_rules! impl_from_uint {
    ($($i:ty),*) => {
        $(
            impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> From<$i> for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
            where
                [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
                [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:,
                [(); EXP_BASE - 2]:
            {
                fn from(n: $i) -> Self
                {
                    Self::from_uint(n)
                }
            }
        )*
    };
}

impl_from_uint!(u8, u16, u32, usize, u64, u128);

/*#[cfg(feature = "ethnum")]
impl_from_uint!(ethnum::U256);*/

macro_rules! impl_from_int {
    ($($i:ty),*) => {
        $(
            impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> From<$i> for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
            where
                [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
                [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:,
                [(); EXP_BASE - 2]:
            {
                fn from(n: $i) -> Self
                {
                    Self::from_int(n)
                }
            }
        )*
    };
}

impl_from_int!(i8, i16, i32, isize, i64, i128);

/*#[cfg(feature = "ethnum")]
impl_from_int!(ethnum::I256);*/

impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> From<f32> for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:,
    [(); EXP_BASE - 2]:
{
    fn from(f: f32) -> Self
    {
        Self::from_fp::<u32, 8, 0, 23, 2>(FpSingle::from_bits(f.to_bits()))
    }
}
impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> From<f64> for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:,
    [(); EXP_BASE - 2]:
{
    fn from(f: f64) -> Self
    {
        Self::from_fp::<u64, 11, 0, 52, 2>(FpDouble::from_bits(f.to_bits()))
    }
}