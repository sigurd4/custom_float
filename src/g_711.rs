use crate::Fp;

/// G.711 8-bit minifloat
pub type FpG711 = Fp<u8, true, 3, 0, 4, 2>;

#[cfg(test)]
pub(crate) macro for_floats {
    ($f:ident, $dir:ident, $expr:expr) => {
        crate::tests::for_floats!(
            FpG711 in g_711 => $f, $dir, $expr
        )
    }
}