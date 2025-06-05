use core::f128::consts::*;

use num_traits::FloatConst;

use crate::{util, Fp, UInt};

impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> FloatConst for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:
{
    #[inline]
    #[doc = "Return Euler’s number."]
    fn E() -> Self
    {
        Self::from(E)
    }

    #[inline]
    #[doc = "Return `1.0 / π`."]
    fn FRAC_1_PI() -> Self
    {
        Self::from(FRAC_1_PI)
    }

    #[inline]
    #[doc = "Return `1.0 / sqrt(2.0)`."]
    fn FRAC_1_SQRT_2() -> Self
    {
        Self::from(FRAC_1_SQRT_2)
    }

    #[inline]
    #[doc = "Return `2.0 / π`."]
    fn FRAC_2_PI() -> Self
    {
        Self::from(FRAC_2_PI)
    }

    #[inline]
    #[doc = "Return `2.0 / sqrt(π)`."]
    fn FRAC_2_SQRT_PI() -> Self
    {
        Self::from(FRAC_2_SQRT_PI)
    }

    #[inline]
    #[doc = "Return `π / 2.0`."]
    fn FRAC_PI_2() -> Self
    {
        Self::from(FRAC_PI_2)
    }

    #[inline]
    #[doc = "Return `π / 3.0`."]
    fn FRAC_PI_3() -> Self
    {
        Self::from(FRAC_PI_3)
    }

    #[inline]
    #[doc = "Return `π / 4.0`."]
    fn FRAC_PI_4() -> Self
    {
        Self::from(FRAC_PI_4)
    }

    #[inline]
    #[doc = "Return `π / 6.0`."]
    fn FRAC_PI_6() -> Self
    {
        Self::from(FRAC_PI_6)
    }

    #[inline]
    #[doc = "Return `π / 8.0`."]
    fn FRAC_PI_8() -> Self
    {
        Self::from(FRAC_PI_8)
    }

    #[inline]
    #[doc = "Return `ln(10.0)`."]
    fn LN_10() -> Self
    {
        Self::from(LN_10)
    }

    #[inline]
    #[doc = "Return `ln(2.0)`."]
    fn LN_2() -> Self
    {
        Self::from(LN_2)
    }

    #[inline]
    #[doc = "Return `log10(e)`."]
    fn LOG10_E() -> Self
    {
        Self::from(LOG10_E)
    }

    #[inline]
    #[doc = "Return `log2(e)`."]
    fn LOG2_E() -> Self
    {
        Self::from(LOG2_E)
    }

    #[inline]
    #[doc = "Return Archimedes’ constant `π`."]
    fn PI() -> Self
    {
        Self::from(PI)
    }

    #[inline]
    #[doc = "Return `sqrt(2.0)`."]
    fn SQRT_2() -> Self
    {
        Self::from(SQRT_2)
    }
}