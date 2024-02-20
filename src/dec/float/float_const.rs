use num_traits::FloatConst;

use crate::{util, Dec, UInt};

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> FloatConst for Dec<U, EXP_SIZE, FRAC_SIZE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 2]:
{
    #[inline]
    #[doc = "Return Euler’s number."]
    fn E() -> Self
    {
        Self::from(f64::E())
    }

    #[inline]
    #[doc = "Return `1.0 / π`."]
    fn FRAC_1_PI() -> Self
    {
        Self::from(f64::FRAC_1_PI())
    }

    #[inline]
    #[doc = "Return `1.0 / sqrt(2.0)`."]
    fn FRAC_1_SQRT_2() -> Self
    {
        Self::from(f64::FRAC_1_SQRT_2())
    }

    #[inline]
    #[doc = "Return `2.0 / π`."]
    fn FRAC_2_PI() -> Self
    {
        Self::from(f64::FRAC_2_PI())
    }

    #[inline]
    #[doc = "Return `2.0 / sqrt(π)`."]
    fn FRAC_2_SQRT_PI() -> Self
    {
        Self::from(f64::FRAC_2_SQRT_PI())
    }

    #[inline]
    #[doc = "Return `π / 2.0`."]
    fn FRAC_PI_2() -> Self
    {
        Self::from(f64::FRAC_PI_2())
    }

    #[inline]
    #[doc = "Return `π / 3.0`."]
    fn FRAC_PI_3() -> Self
    {
        Self::from(f64::FRAC_PI_3())
    }

    #[inline]
    #[doc = "Return `π / 4.0`."]
    fn FRAC_PI_4() -> Self
    {
        Self::from(f64::FRAC_PI_4())
    }

    #[inline]
    #[doc = "Return `π / 6.0`."]
    fn FRAC_PI_6() -> Self
    {
        Self::from(f64::FRAC_PI_6())
    }

    #[inline]
    #[doc = "Return `π / 8.0`."]
    fn FRAC_PI_8() -> Self
    {
        Self::from(f64::FRAC_PI_8())
    }

    #[inline]
    #[doc = "Return `ln(10.0)`."]
    fn LN_10() -> Self
    {
        Self::from(f64::LN_10())
    }

    #[inline]
    #[doc = "Return `ln(2.0)`."]
    fn LN_2() -> Self
    {
        Self::from(f64::LN_2())
    }

    #[inline]
    #[doc = "Return `log10(e)`."]
    fn LOG10_E() -> Self
    {
        Self::from(f64::LOG10_E())
    }

    #[inline]
    #[doc = "Return `log2(e)`."]
    fn LOG2_E() -> Self
    {
        Self::from(f64::LOG2_E())
    }

    #[inline]
    #[doc = "Return Archimedes’ constant `π`."]
    fn PI() -> Self
    {
        Self::from(f64::PI())
    }

    #[inline]
    #[doc = "Return `sqrt(2.0)`."]
    fn SQRT_2() -> Self
    {
        Self::from(f64::SQRT_2())
    }
}