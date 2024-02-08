use num_traits::FloatConst;

use crate::{bitsize_of, Fp, UInt};

impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> FloatConst for Fp<U, EXP_SIZE, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 1]:
{
    #[doc = "Return Euler’s number."]
    fn E() -> Self
    {
        Self::from(f64::E())
    }

    #[doc = "Return `1.0 / π`."]
    fn FRAC_1_PI() -> Self
    {
        Self::from(f64::FRAC_1_PI())
    }

    #[doc = "Return `1.0 / sqrt(2.0)`."]
    fn FRAC_1_SQRT_2() -> Self
    {
        Self::from(f64::FRAC_1_SQRT_2())
    }

    #[doc = "Return `2.0 / π`."]
    fn FRAC_2_PI() -> Self
    {
        Self::from(f64::FRAC_2_PI())
    }

    #[doc = "Return `2.0 / sqrt(π)`."]
    fn FRAC_2_SQRT_PI() -> Self
    {
        Self::from(f64::FRAC_2_SQRT_PI())
    }

    #[doc = "Return `π / 2.0`."]
    fn FRAC_PI_2() -> Self
    {
        Self::from(f64::FRAC_PI_2())
    }

    #[doc = "Return `π / 3.0`."]
    fn FRAC_PI_3() -> Self
    {
        Self::from(f64::FRAC_PI_3())
    }

    #[doc = "Return `π / 4.0`."]
    fn FRAC_PI_4() -> Self
    {
        Self::from(f64::FRAC_PI_4())
    }

    #[doc = "Return `π / 6.0`."]
    fn FRAC_PI_6() -> Self
    {
        Self::from(f64::FRAC_PI_6())
    }

    #[doc = "Return `π / 8.0`."]
    fn FRAC_PI_8() -> Self
    {
        Self::from(f64::FRAC_PI_8())
    }

    #[doc = "Return `ln(10.0)`."]
    fn LN_10() -> Self
    {
        Self::from(f64::LN_10())
    }

    #[doc = "Return `ln(2.0)`."]
    fn LN_2() -> Self
    {
        Self::from(f64::LN_2())
    }

    #[doc = "Return `log10(e)`."]
    fn LOG10_E() -> Self
    {
        Self::from(f64::LOG10_E())
    }

    #[doc = "Return `log2(e)`."]
    fn LOG2_E() -> Self
    {
        Self::from(f64::LOG2_E())
    }

    #[doc = "Return Archimedes’ constant `π`."]
    fn PI() -> Self
    {
        Self::from(f64::PI())
    }

    #[doc = "Return `sqrt(2.0)`."]
    fn SQRT_2() -> Self
    {
        Self::from(f64::SQRT_2())
    }
}