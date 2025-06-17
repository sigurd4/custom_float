use crate::Fp;

/// Khronos 10-bit unsigned mini-float
pub type KhronosFp10 = Fp<u16, false, 5, 0, 5, 2>;

/// Khronos 11-bit unsigned mini-float
pub type KhronosFp11 = Fp<u16, false, 5, 0, 6, 2>;

#[cfg(test)]
pub(crate) macro for_floats {
    ($f:ident, $dir:ident, $expr:expr) => {
        crate::tests::for_floats!(
            KhronosFp10 in khronos,
            KhronosFp11 in khronos => $f, $dir, $expr
        )
    }
}