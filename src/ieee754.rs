//IEEE754

use crate::Fp;

/// Half-precision floating-point format, known as binary16, FP16 or float16
pub type FpHalf = Fp<u16, true, 5, 0, 10, 2>;

/// Single-precision floating-point format, known as binary32, FP32 or float32
pub type FpSingle = Fp<u32, true, 8, 0, 23, 2>;

/// Double-precision floating-point format, known as binary64, FP64 or float64
pub type FpDouble = Fp<u64, true, 11, 0, 52, 2>;

/// Quadruple-precision floating-point format, known as binary128, FP128 or float128
pub type FpQuadruple = Fp<u128, true, 15, 0, 112, 2>;

// /// Octuple-precision floating-point format, known as binary256, FP256 or float256
//pub type FpOctuple = Fp<u256, true, 19, 0, 236, 2>;

/// Decimal single-precision floating-point format, known as decimal32
pub type DecSingle = Fp<u32, true, 11, 7, 13, 10>;

/// Decimal double-precision floating-point format, known as decimal64
pub type DecDouble = Fp<u64, true, 13, 16, 34, 10>;

/// Decimal quadruple-precision floating-point format, known as decimal128
pub type DecQuadruple = Fp<u128, true, 17, 34, 76, 10>;

#[cfg(test)]
pub(crate) macro for_floats {
    ($f:ident, $dir:ident, $expr:expr) => {
        crate::tests::for_floats!(
            FpHalf in ieee754,
            FpSingle in ieee754,
            FpDouble in ieee754,
            FpQuadruple in ieee754,
            //FpOctuple in ieee754,
            DecSingle in ieee754,
            DecDouble in ieee754,
            DecQuadruple in ieee754 => $f, $dir, $expr
        )
    }
}