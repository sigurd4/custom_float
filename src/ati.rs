use crate::Fp;

/// Radeon R300 / R420 FP24
pub type Fp24 = Fp<u32, true, 7, 0, 16, 2>;

#[cfg(test)]
pub(crate) macro for_floats {
    ($f:ident, $dir:ident, $expr:expr) => {
        crate::tests::for_floats!(
            Fp24 in ati => $f, $dir, $expr
        )
    }
}