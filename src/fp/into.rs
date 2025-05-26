use crate::{ieee754::{FpDouble, FpHalf, FpQuadruple, FpSingle}, util, Fp, UInt};

macro_rules! impl_into_float {
    ($($f:ty: $fp:ty),*) => {
        $(
            impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> From<Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>> for $f
            where
                [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
                [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - 0 - FRAC_SIZE]:,
                [(); EXP_BASE - 2]:
            {
                fn from(value: Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>) -> Self
                {
                    <$f>::from_bits(<$fp>::from_fp(value).to_bits())
                }
            }
        )*
    };
}
impl_into_float!(f16: FpHalf, f32: FpSingle, f64: FpDouble, f128: FpQuadruple);