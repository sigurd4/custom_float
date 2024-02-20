use num_traits::AsPrimitive;

use crate::{util, Fp, UInt};

macro_rules! impl_as_primitive_uint {
    ($($i:ty),*) => {
        $(
            impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> AsPrimitive<$i> for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
            where
                U: 'static,
                [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
                [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:,
                [(); EXP_BASE - 2]:
            {
                #[inline]
                fn as_(self) -> $i
                {
                    self.to_uint_wrapping()
                }
            }
        )*
    };
}

impl_as_primitive_uint!(u8, u16, u32, usize, u64, u128/*, U256*/);

macro_rules! impl_as_primitive_int {
    ($($i:ty),*) => {
        $(
            impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> AsPrimitive<$i> for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
            where
                U: 'static,
                [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
                [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:,
                [(); EXP_BASE - 2]:
            {
                #[inline]
                fn as_(self) -> $i
                {
                    self.to_int_wrapping()
                }
            }
        )*
    };
}

impl_as_primitive_int!(i8, i16, i32, isize, i64, i128/*, I256*/);

impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> AsPrimitive<f32> for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: 'static,
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:,
    [(); EXP_BASE - 2]:
{
    #[inline]
    fn as_(self) -> f32
    {
        Into::<f32>::into(self)
    }
}

impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> AsPrimitive<f64> for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: 'static,
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:,
    [(); EXP_BASE - 2]:
{
    #[inline]
    fn as_(self) -> f64
    {
        Into::<f64>::into(self)
    }
}