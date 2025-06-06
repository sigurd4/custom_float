use crate::{FpRepr, ieee754::{FpDouble, FpHalf, FpQuadruple, FpSingle}, Fp};

macro_rules! impl_from_uint {
    ($($i:ty),*) => {
        $(
            impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> From<$i> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
            where
                U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
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
            impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> From<$i> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
            where
                U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
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

macro_rules! impl_from_float {
    ($($f:ty: $b:ty, $fp:ty),*) => {
        $(
            impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> From<$f> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
            where
                U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
            {
                fn from(f: $f) -> Self
                {
                    Self::from_fp::<$b, _, _, _, _, 2>(<$fp>::from_bits(f.to_bits()))
                }
            }
        )*
    };
}
impl_from_float!(f16: u16, FpHalf, f32: u32, FpSingle, f64: u64, FpDouble, f128: u128, FpQuadruple);