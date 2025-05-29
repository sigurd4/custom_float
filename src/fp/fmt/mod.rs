moddef::moddef!(
    mod {
        binary,
        debug,
        display,
        lower_exp
    }
);

macro_rules! suited {
    ($t:ty) => {
        <$t as From<Self>>::from(Self::min_value()) != 0.0
    };
}
macro_rules! do_fmt {
    ($value:expr, $f:expr => $trait:ident, $t:ty) => {
        {
            let x = <$t as From<Self>>::from(*$value);
            if $value.is_finite() != x.is_finite()
            {
                write!($f, "~")?
            }
            return <$t as $trait>::fmt(&x, $f)
        }
    };
}
macro_rules! do_maybe_fmt {
    ($value:expr, $f:expr => $trait:ident, $t:ty) => {
        {
            let x = <$t as From<Self>>::from(*$value);
            if $value.is_finite() == x.is_finite()
            {
                return <$t as $trait>::fmt(&x, $f)
            }
        }
    };
}

macro_rules! fmt_as {
    ($value:expr, $f:expr => $trait:ident => $or:ty) => {
        {
            if super::suited!(f16)
            {
                super::do_maybe_fmt!($value, $f => $trait, f16)
            }
            if super::suited!(f32)
            {
                super::do_maybe_fmt!($value, $f => $trait, f32)
            }
            if super::suited!(f64)
            {
                super::do_maybe_fmt!($value, $f => $trait, f64)
            }
            super::do_fmt!($value, $f => $trait, $or)
        }
    };
}

use fmt_as as fmt_as;
use do_fmt as do_fmt;
use do_maybe_fmt as do_maybe_fmt;
use suited as suited;