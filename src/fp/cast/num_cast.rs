use num_traits::{ToPrimitive, NumCast};

use crate::{util, Fp, UInt};

trait NumCastSpec1: ToPrimitive
{
    fn __cast<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize>(self) -> Option<Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>>
    where
        [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
        [(); EXP_BASE - 2]:;
}

impl<T> NumCastSpec1 for T
where
    T: ToPrimitive
{
    default fn __cast<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize>(self) -> Option<Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>>
    where
        [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
        [(); EXP_BASE - 2]:
    {
        self.to_f64().map(|f| f.into())
    }
}
impl<T: UInt, const S: bool, const E: usize, const I: usize, const F: usize, const B: usize> NumCastSpec1 for Fp<T, S, E, I, F, B>
where
    [(); util::bitsize_of::<T>() - S as usize - E - I - F]:,
    [(); B - 2]:
{
    fn __cast<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize>(self) -> Option<Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>>
    where
        [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
        [(); EXP_BASE - 2]:
    {
        Some(Fp::from_fp(self))
    }
}

trait NumCastSpec<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize>: ToPrimitive
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:
{
    fn _cast(self) -> Option<Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>>;
}
impl<T, U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> NumCastSpec<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE> for T
where
    T: ToPrimitive,
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:
{
    default fn _cast(self) -> Option<Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>>
    {
        NumCastSpec1::__cast(self)
    }
}
impl<T, U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> NumCastSpec<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE> for T
where
    T: ToPrimitive + Into<Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>>,
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:
{
    fn _cast(self) -> Option<Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>>
    {
        Some(self.into())
    }
}

impl<U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> NumCast for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    [(); util::bitsize_of::<U>() - SIGN_BIT as usize - EXP_SIZE - INT_SIZE - FRAC_SIZE]:,
    [(); EXP_BASE - 2]:
{
    #[inline]
    fn from<T: ToPrimitive>(n: T) -> Option<Self>
    {
        n._cast()
    }
}