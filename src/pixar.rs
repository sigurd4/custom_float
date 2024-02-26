use crate::Fp;

/// PXR24 compressed floating-point format
pub type FpPXR24 = Fp<u32, true, 8, 0, 15, 2>;