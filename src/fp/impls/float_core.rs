use num_traits::{float::FloatCore, Float};

use crate::{bitsize_of, Fp, UInt};


impl<U: UInt, const EXP_SIZE: usize, const INT_BIT: bool, const FRAC_SIZE: usize> FloatCore for Fp<U, EXP_SIZE, INT_BIT, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - INT_BIT as usize - FRAC_SIZE - 1]:,
    [(); bitsize_of::<U>() - EXP_SIZE - false as usize - FRAC_SIZE - 1]:
{
    fn infinity() -> Self
    {
        <Self as Float>::infinity()
    }

    fn neg_infinity() -> Self
    {
        <Self as Float>::neg_infinity()
    }

    fn nan() -> Self
    {
        <Self as Float>::nan()
    }

    fn neg_zero() -> Self
    {
        <Self as Float>::neg_zero()
    }

    fn min_value() -> Self
    {
        <Self as Float>::min_value()
    }

    fn min_positive_value() -> Self
    {
        <Self as Float>::min_positive_value()
    }

    fn epsilon() -> Self
    {
        <Self as Float>::epsilon()
    }

    fn max_value() -> Self
    {
        <Self as Float>::max_value()
    }

    fn classify(self) -> std::num::FpCategory
    {
        <Self as Float>::classify(self)
    }

    fn to_degrees(self) -> Self
    {
        <Self as Float>::to_degrees(self)
    }

    fn to_radians(self) -> Self
    {
        <Self as Float>::to_radians(self)
    }

    fn integer_decode(self) -> (u64, i16, i8)
    {
        <Self as Float>::integer_decode(self)
    }
}