use num_traits::{float::FloatCore, Float, NumCast};

use crate::{bitsize_of, Fp, UInt};


impl<U: UInt, const EXP_SIZE: usize, const INT_BIT: bool, const FRAC_SIZE: usize> FloatCore for Fp<U, EXP_SIZE, INT_BIT, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - INT_BIT as usize - FRAC_SIZE - 1]:,
    [(); bitsize_of::<U>() - EXP_SIZE - false as usize - FRAC_SIZE - 1]:
{
    fn infinity() -> Self
    {
        Self::infinity()
    }

    fn neg_infinity() -> Self
    {
        Self::neg_infinity()
    }

    fn nan() -> Self
    {
        Self::nan()
    }

    fn neg_zero() -> Self
    {
        Self::neg_zero()
    }

    fn min_value() -> Self
    {
        Self::min_value()
    }

    fn min_positive_value() -> Self
    {
        Self::min_positive_value()
    }

    fn epsilon() -> Self
    {
        Self::epsilon()
    }

    fn max_value() -> Self
    {
        Self::max_value()
    }

    fn classify(self) -> std::num::FpCategory
    {
        self.classify()
    }

    fn to_degrees(self) -> Self
    {
        self.to_degrees()
    }

    fn to_radians(self) -> Self
    {
        self.to_radians()
    }

    fn integer_decode(self) -> (u64, i16, i8)
    {
        let s = self.sign_bit();
        let e = self.exp_bits();
        let mut f = self.frac_bits();
        let bias = Self::exp_bias();

        if INT_BIT
        {
            f = f + (self.int_bit() << FRAC_SIZE)
        }

        let s = NumCast::from(s).unwrap();
        let e = NumCast::from(<i32 as NumCast>::from(e).unwrap() - <i32 as NumCast>::from(bias).unwrap()).unwrap();
        let f = NumCast::from(f).unwrap();

        (f, e, s)
    }
}