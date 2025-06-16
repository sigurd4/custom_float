use crate::Fp;

/// PXR24 compressed floating-point format
pub type FpPXR24 = Fp<u32, true, 8, 0, 15, 2>;

#[cfg(test)]
pub(crate) macro for_floats {
    ($expr:expr) => {
        crate::tests::for_floats!(
            FpPXR24 in pixar => $expr
        )
    }
}