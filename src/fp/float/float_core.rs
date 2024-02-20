use num_traits::{float::FloatCore, NumCast};

use crate::{util, Fp, UInt};


impl<U: UInt, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> FloatCore for Fp<U, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - EXP_SIZE - INT_SIZE - FRAC_SIZE - 1]:,
    [(); util::bitsize_of::<U>() - EXP_SIZE - 0 - FRAC_SIZE - 1]:,
    [(); EXP_BASE - 2]:
{
    #[inline]
    fn infinity() -> Self
    {
        Self::infinity()
    }

    #[inline]
    fn neg_infinity() -> Self
    {
        Self::neg_infinity()
    }

    #[inline]
    fn nan() -> Self
    {
        Self::nan()
    }

    #[inline]
    fn neg_zero() -> Self
    {
        Self::neg_zero()
    }

    #[inline]
    fn min_value() -> Self
    {
        Self::min_value()
    }

    #[inline]
    fn min_positive_value() -> Self
    {
        Self::min_positive_value()
    }

    #[inline]
    fn epsilon() -> Self
    {
        Self::epsilon()
    }

    #[inline]
    fn max_value() -> Self
    {
        Self::max_value()
    }

    #[inline]
    fn classify(self) -> std::num::FpCategory
    {
        self.classify()
    }

    #[inline]
    fn to_degrees(self) -> Self
    {
        self.to_degrees()
    }

    #[inline]
    fn to_radians(self) -> Self
    {
        self.to_radians()
    }

    #[inline]
    fn integer_decode(self) -> (u64, i16, i8)
    {
        let s = self.sign_bit();
        let e = self.exp_bits();
        let mut f = self.frac_bits();
        let bias = Self::exp_bias();

        if INT_SIZE > 0
        {
            f = f + (self.int_bits() << FRAC_SIZE)
        }

        let s = 1 - <i8 as NumCast>::from(s).unwrap()*2;
        let e = NumCast::from(<i32 as NumCast>::from(e).unwrap() - <i32 as NumCast>::from(bias).unwrap()).unwrap();
        let f = NumCast::from(f).unwrap();

        (f, e, s)
    }
}