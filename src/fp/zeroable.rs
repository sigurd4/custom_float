use bytemuck::{Zeroable, ZeroableInOption};

use crate::{Fp, FpRepr};

unsafe impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> ZeroableInOption for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    Self: Zeroable,
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE> + 'static,
{

}

#[cfg(test)]
mod test
{
    use bytemuck::Zeroable;

    use crate::tests::for_all_floats;

    #[test]
    fn zeroable()
    {
        for_all_floats!(F, _DIR, {
            let x = F::zeroed();
            assert!(x.is_zero())
        })
    }

    #[test]
    fn zeroable_in_option()
    {
        for_all_floats!(F, _DIR, {
            let x = Option::<F>::zeroed();
            assert_eq!(x, None)
        })
    }
}