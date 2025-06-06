use crate::{FpRepr, ieee754::{FpDouble, FpHalf, FpQuadruple, FpSingle}, Fp};

macro_rules! impl_into_float {
    ($($f:ty: $fp:ty),*) => {
        $(
            impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> From<Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>> for $f
            where
                U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
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