use crate::Fp;

/// x86 extended precision floating-point format
pub type Fp80 = Fp<u128, true, 15, 1, 63, 2>;

#[cfg(test)]
pub(crate) macro for_floats {
    ($f:ident, $dir:ident, $expr:expr) => {
        crate::tests::for_floats!(
            Fp80 in intel => $f, $dir, $expr
        )
    }
}