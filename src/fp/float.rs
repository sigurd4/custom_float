use std::{cmp::Ordering, num::FpCategory};

use array_math::{ArrayMath, ArrayOps};
use num_traits::{Float, Inv, FloatConst, NumCast, One, Zero};

use super::{bitsize_of, Fp, UInt};


impl<U: UInt, const EXP_SIZE: usize, const FRAC_SIZE: usize> Float for Fp<U, EXP_SIZE, FRAC_SIZE>
where
    [(); bitsize_of::<U>() - EXP_SIZE - FRAC_SIZE - 1]:
{
    fn nan() -> Self
    {
        Self::nan()
    }

    fn infinity() -> Self
    {
        Self::infinity()
    }

    fn neg_infinity() -> Self
    {
        -Self::infinity()
    }

    fn neg_zero() -> Self
    {
        -Self::zero()
    }

    fn min_value() -> Self
    {
        -Self::max_value()
    }

    fn min_positive_value() -> Self
    {
        Self::from_bits(U::one())
    }

    fn max_value() -> Self
    {
        Self::from_bits(((U::max_value() >> bitsize_of::<U>() - EXP_SIZE) << Self::EXP_POS) - U::one())
    }

    fn is_nan(self) -> bool
    {
        !self.is_finite() && !self.frac_bits().is_zero()
    }

    fn is_infinite(self) -> bool
    {
        !self.is_finite() && self.frac_bits().is_zero()
    }

    fn is_finite(self) -> bool
    {
        self.exp_bits() != (U::max_value() >> bitsize_of::<U>() - EXP_SIZE)
    }

    fn is_normal(self) -> bool
    {
        let e = self.exp_bits();
        e != U::zero() && e != (U::max_value() >> bitsize_of::<U>() - EXP_SIZE)
    }

    fn classify(self) -> FpCategory
    {
        let e = self.exp_bits();
        if e == U::zero()
        {
            if self.is_zero()
            {
                return FpCategory::Zero
            }
            return FpCategory::Subnormal
        }
        if e == (U::max_value() >> bitsize_of::<U>() - EXP_SIZE)
        {
            if self.frac_bits().is_zero()
            {
                return FpCategory::Infinite
            }
            return FpCategory::Nan
        }
        FpCategory::Normal
    }

    fn floor(self) -> Self
    {
        if self.is_zero()
        {
            return Self::zero()
        }
        let s = !self.sign_bit().is_zero();
        let e = self.exp_bits();
        let mut f = self.frac_bits();
        let bias = Self::exp_bias();
        if e < bias
        {
            return if s {-Self::one()} else {Self::zero()}
        }
        let mask = match ((U::one() << FRAC_SIZE) - U::one()).checked_shr(<u32 as NumCast>::from(e - bias).unwrap())
        {
            Some(mask) => mask,
            None => U::zero()
        };
        let i = s && !(f & mask).is_zero();
        f = f & (!mask);
        let s_bit = if s {U::one() << Self::SIGN_POS} else {U::zero()};
        let n = Self::from_bits(s_bit + (e << Self::EXP_POS) + f);
        
        if i
        {
            return n - Self::one()
        }
        n
    }

    fn ceil(self) -> Self
    {
        if self.is_zero()
        {
            return Self::zero()
        }
        let s = !self.sign_bit().is_zero();
        let e = self.exp_bits();
        let mut f = self.frac_bits();
        let bias = Self::exp_bias();
        if e < bias
        {
            return if s {Self::zero()} else {Self::one()}
        }
        let mask = match ((U::one() << FRAC_SIZE) - U::one()).checked_shr(<u32 as NumCast>::from(e - bias).unwrap())
        {
            Some(mask) => mask,
            None => U::zero()
        };
        let i = !s && !(f & mask).is_zero();
        f = f & (!mask);
        let s_bit = if s {U::one() << Self::SIGN_POS} else {U::zero()};
        let n = Self::from_bits(s_bit + (e << Self::EXP_POS) + f);
        
        if i
        {
            return n + Self::one()
        }
        n
    }

    fn round(self) -> Self
    {
        let c = self.ceil();
        let f = self.floor();
        if if self >= Self::zero() {c - self <= self - f} else {c - self < self - f}
        {
            c
        }
        else
        {
            f
        }
    }

    fn trunc(self) -> Self
    {
        let s = !self.sign_bit().is_zero();
        if self.is_zero()
        {
            return if s {-Self::zero()} else {Self::zero()}
        }
        let e = self.exp_bits();
        let mut f = self.frac_bits();
        let bias = Self::exp_bias();
        if e < bias
        {
            return if s {-Self::zero()} else {Self::zero()}
        }
        let mask = match ((U::one() << FRAC_SIZE) - U::one()).checked_shr(<u32 as NumCast>::from(e - bias).unwrap())
        {
            Some(mask) => mask,
            None => U::zero()
        };
        f = f & (!mask);
        let s_bit = if s {U::one() << Self::SIGN_POS} else {U::zero()};
        return Self::from_bits(s_bit + (e << Self::EXP_POS) + f);
    }

    fn fract(self) -> Self
    {
        self - self.trunc()
    }

    fn abs(self) -> Self
    {
        return Self::from_bits(self.to_bits() & (U::max_value() >> bitsize_of::<U>() - Self::SIGN_POS));
    }

    fn signum(self) -> Self
    {
        if self.is_nan()
        {
            return self
        }
        let s = !self.sign_bit().is_zero();
        if s {-Self::one()} else {Self::one()}
    }

    fn is_sign_positive(self) -> bool
    {
        self.sign_bit().is_zero()
    }

    fn is_sign_negative(self) -> bool
    {
        !self.sign_bit().is_zero()
    }

    fn mul_add(self, a: Self, b: Self) -> Self
    {
        (self*a) + b
    }

    fn recip(self) -> Self
    {
        self.inv()
    }

    fn powi(self, mut n: i32) -> Self
    {
        let mut r = Self::one();
        
        let mut x = if n < 0
        {
            n = -n;
            Self::one() / self
        }
        else
        {
            self
        };
        
        if (n & 1) == 1
        {
            r = x;
        }

        loop
        {
            n >>= 1;
            if 0 == n
            {
                break;
            }
            x = x*x;
            if (n & 1) == 1
            {
                r = r*x;
            }
        }

        return r;
    }

    fn powf(self, n: Self) -> Self
    {
        let xabs = self.abs();
        let nabs = n.abs();

        let edge_x = {
            let e = (U::one() << EXP_SIZE) - U::one();
            Self::from_bits(e << Self::EXP_POS)
        };
        let edge_n = {
            let e = Self::exp_bias() + U::from(FRAC_SIZE).unwrap();
            Self::from_bits(e << Self::EXP_POS)
        };
        
        if n.is_zero()
        {
            return Self::one()
        }
        if self.is_one() && !n.is_finite()
        {
            return Self::one()
        }
        if self.is_nan() || n.is_nan()
        {
            return self + n
        }
        if n == Self::one()
        {
            return self
        }
        
        // if n is an odd integer
        let n_d = (nabs + Self::one())/<Self as From<_>>::from(2.0);
        let noi = n_d.trunc() == n_d;

        if xabs.is_zero()
        {
            if noi || n.is_infinite()
            {
                if n < Self::zero()
                {
                    return Self::one()/self
                }
                return self
            }
            if n > Self::zero()
            {
                return Self::zero()
            }
            return Self::one()/xabs
        }

        if (xabs.to_bits() >= edge_x.to_bits()) || (nabs.to_bits() >= edge_n.to_bits())
        {
            if nabs.is_infinite()
            {
                if self == -Self::one()
                {
                    return Self::one()
                }
                if xabs > Self::one()
                {
                    if n > Self::zero()
                    {
                        return n
                    }
                    return Self::zero()
                }
                if n > Self::zero()
                {
                    return Self::zero()
                }
                return nabs
            }

            if self == Self::infinity()
            {
                if n < Self::zero()
                {
                    return Self::zero()
                }
                return self
            }
            if self > Self::neg_infinity()
            {
                if !noi
                {
                    if self < Self::zero() || n > Self::zero()
                    {
                        return self.sqrt()
                    }
                    return (Self::one()/self).sqrt()
                }
            
                // clamp  -0x1.0p31 < y < 0x1.0p31
                let nabs = nabs.min(NumCast::from(U::max_value()).unwrap());
                let s = !n.sign_bit().is_zero();
                let mut i = <u128 as NumCast>::from(nabs).unwrap();
                let mut r = Self::one();
                
                let mut x = if s
                {
                    Self::one() / self
                }
                else
                {
                    self
                };
                
                if (i & 1) == 1
                {
                    r = x;
                }

                loop
                {
                    i >>= 1;
                    if 0 == i
                    {
                        break;
                    }
                    x = x*x;
                    if (i & 1) == 1
                    {
                        r = r*x;
                    }
                }

                return r;
            }

            if noi
            {
                if n > Self::zero()
                {
                    return self
                }
                return -Self::zero()
            }
            
            if n > Self::zero()
            {
                return -self
            }
            return Self::zero()
        }
        
        let ni = n.trunc() == n;

        if !ni && self < Self::zero()
        {
            return Self::nan()
        }
        
        //special case for sqrts
        if <Self as From<_>>::from(0.5) == nabs
        {
            if self < Self::zero() || n > Self::zero()
            {
                return self.sqrt()
            }
            return (Self::one()/self).sqrt()
        }

        let xabs_log2 = xabs.log2();

        let n_x_abs_log2 = n*xabs_log2;
        if n_x_abs_log2.is_nan()
        {
            println!("NaN");
        }
        if n_x_abs_log2 >= NumCast::from(U::one() << (EXP_SIZE - 1)).unwrap()
        {
            // ???
            return Self::from_bits(((U::one() << EXP_SIZE) - U::one()) << Self::EXP_POS)*n_x_abs_log2;
        }
        if n_x_abs_log2 <= -<Self as NumCast>::from((U::one() << (EXP_SIZE - 1)) + U::from(FRAC_SIZE - 1).unwrap()).unwrap()
        {
            // ???
            return -Self::from_bits(U::one() << (Self::EXP_POS - EXP_SIZE/2))*n_x_abs_log2;
        }
        
        let e = n_x_abs_log2.floor();
        let f = n_x_abs_log2 - e;
        
        const C0: f64 =  1.0 + 0.278626872016317130037181614004e-10;
        const C1: f64 = 0.693147176943623740308984004029708;
        const C2: f64 = 0.240226505817268621584559118975830;
        const C3: f64 = 0.555041568519883074165425891257052e-1;
        const C4_OVER_C8: f64 = 0.961813690023115610862381719985771e-2 / C8;
        const C5_OVER_C8: f64 = 0.133318252930790403741964203236548e-2 / C8;
        const C6_OVER_C8: f64 = 0.154016177542147239746127455226575e-3 / C8;
        const C7_OVER_C8: f64 = 0.154832722143258821052933667742417e-4 / C8;
        const C8: f64 = 0.134107709538786543922336536865157e-5;

        let mut z = Self::one();
        if !f.is_zero()
        {
            let ff = f * f;
            let mut s7 = <Self as From<_>>::from(C7_OVER_C8) * f;
            let mut s3 = <Self as From<_>>::from(C3) * f;
            let mut s5 = <Self as From<_>>::from(C5_OVER_C8) * f;
            let mut s1 = <Self as From<_>>::from(C1) * f;
            let ffff = ff * ff;
            s7 += <Self as From<_>>::from(C6_OVER_C8);
            s3 += <Self as From<_>>::from(C2);
            s5 += <Self as From<_>>::from(C4_OVER_C8);
            s1 += <Self as From<_>>::from(C0);
            s7 *= ff;
            s3 *= ff;
            s5 += ffff;
            let c8ffff = ffff * <Self as From<_>>::from(C8);
            s7 += s5;
            s3 += s1;
            s7 *= c8ffff;
            z = s3 + s7;
        }

        let e = e + <Self as NumCast>::from(Self::exp_bias()).unwrap();
        if e > <Self as NumCast>::from((U::one() << EXP_SIZE) - U::one()).unwrap()
        {
            return Self::infinity()*z
        }
        if e < Self::zero()
        {
            return Self::zero()*z
        }
        if e.is_nan()
        {
            return Self::nan()
        }
        let e = U::from(e).unwrap();
        let y = Self::from_bits(e << Self::EXP_POS);

        y*z
    }

    fn sqrt(self) -> Self
    {
        if self.is_nan() || self.is_sign_negative()
        {
            return Self::nan()
        }

        let y = Self::from_bits(
            ((Self::exp_bias() + U::one()) * (U::one() << (FRAC_SIZE - 1)))
            + (self.to_bits() >> 1usize)
            - (U::one() << (FRAC_SIZE - 1))
        );

        const NEWTON: usize = 4;
        let half = <Self as From<_>>::from(0.5);
        let mut y = y;
        for _ in 0..NEWTON
        {
            y = half*(y + self/y);
        }
        y
    }

    fn exp(self) -> Self
    {
        (self/Self::LN_2()).exp2()
    }

    fn exp2(self) -> Self
    {
        if self.is_nan()
        {
            return self
        }
        if self.is_infinite()
        {
            if self.is_sign_negative()
            {
                return Self::zero()
            }
            return self
        }
        let neg = self.is_sign_negative();
        let x = self.abs().max(-<Self as NumCast>::from(Self::exp_bias() - U::one()).unwrap());
        let w = if neg {x.floor()} else {x.floor()};
        let z = x - w;
        let y = <Self as NumCast>::from(U::one() << FRAC_SIZE).unwrap()*(x
            + <Self as From<_>>::from(121.2740575)
            + <Self as From<_>>::from(27.7280233)/(<Self as From<_>>::from(4.84252568) - z)
            - <Self as From<_>>::from(1.49012907)*z
        );
        if y.is_nan()
        {
            return y
        }
        if y >= <Self as NumCast>::from(U::one() << Self::BIT_SIZE - 1).unwrap()
        {
            if neg
            {
                return Self::zero()
            }
            return Self::infinity()
        }
        if y <= Self::one()
        {
            if neg
            {
                return Self::infinity()
            }
            return Self::zero();
        }
        let y = Self::from_bits(U::from(y).unwrap());
        if neg
        {
            return y.inv();
        }
        y
    }

    fn ln(self) -> Self
    {
        self.log2()*Self::LN_2()
    }

    fn log(self, base: Self) -> Self
    {
        self.ln()/base.ln()
    }

    fn log2(self) -> Self
    {
        if self.is_nan()
        {
            return self
        }
        if self < Self::zero()
        {
            return Self::nan()
        }
        if self.is_zero()
        {
            return Self::neg_infinity()
        }
        if self.is_infinite()
        {
            return Self::infinity()
        }
        let b = if !self.exp_bits().is_zero()
        {
            U::zero()
        }
        else
        {
            U::from(self.frac_bits().leading_zeros() as usize - (bitsize_of::<U>() - FRAC_SIZE)).unwrap()
        };
        let mut y = <Self as NumCast>::from(self.exp_bits()).unwrap()
            - <Self as NumCast>::from(Self::exp_bias() + U::one()).unwrap()
            - <Self as NumCast>::from(b).unwrap();
        let mut u = self.to_bits();
        u = u & !(((U::one() << EXP_SIZE) - U::one()) << Self::EXP_POS);
        u = u + (Self::exp_bias() << Self::EXP_POS);
        let u = Self::from_bits(u);
        y += (<Self as From<_>>::from(-0.34484843) * u + <Self as From<_>>::from(2.02466578)) * u  - <Self as From<_>>::from(0.67487759); 
        return y;
    }

    fn log10(self) -> Self
    {
        self.ln()/Self::LN_10()
    }

    fn max(self, other: Self) -> Self
    {
        match self.partial_cmp(&other)
        {
            Some(cmp) => match cmp
            {
                Ordering::Less => other,
                Ordering::Equal => self,
                Ordering::Greater => self,
            }
            None => self + other
        }
    }

    fn min(self, other: Self) -> Self
    {
        match self.partial_cmp(&other)
        {
            Some(cmp) => match cmp
            {
                Ordering::Less => self,
                Ordering::Equal => self,
                Ordering::Greater => other,
            }
            None => self + other
        }
    }

    fn abs_sub(self, other: Self) -> Self
    {
        (self - other).abs()
    }

    fn cbrt(self) -> Self
    {
        self.powf(Self::one()/<Self as From<_>>::from(3.0))
    }

    fn hypot(self, other: Self) -> Self
    {
        (self*self + other*other).sqrt()
    }

    fn sin(self) -> Self
    {
        const N: usize = 6;
        const C: [f64; N] = [
            1.276278962,
            -0.285261569,
            0.009118016,
            -0.000136587,
            0.000001185,
            -0.000000007
        ];
        let t: [[Self; N]; N] = ArrayOps::fill(
            |n| <[_; N]>::chebyshev_polynomial(1, n).unwrap()
        );
        let p: [Self; N] = t.zip(C)
            .map2(|(t, c)| t.map2(|tn| <Self as From<_>>::from(c)*tn))
            .reduce(|a, b| a.zip(b).map2(|(a, b)| a + b))
            .unwrap_or_default();

        let mut w = self*Self::FRAC_2_PI();
        let mut i = 0;
        while i < 4
        {
            w -= Self::one();
            if i % 2 == 0 && w < Self::zero()
            {
                w = -w;
            }
            w %= <Self as From<_>>::from(4.0);
            i += 1;
        }
        let two = <Self as From<_>>::from(2.0);
        let w = if w > Self::one() {two - w} else if w < -Self::one() {-two - w} else {w};

        let z = two*w*w - Self::one();

        p.polynomial(z)*w
    }

    fn cos(self) -> Self
    {
        const N: usize = 6;
        const C: [f64; N] = [
            1.276278962,
            -0.285261569,
            0.009118016,
            -0.000136587,
            0.000001185,
            -0.000000007
        ];
        let t: [[Self; N]; N] = ArrayOps::fill(
            |n| <[_; N]>::chebyshev_polynomial(1, n).unwrap()
        );
        let p: [Self; N] = t.zip(C)
            .map2(|(t, c)| t.map2(|tn| <Self as From<_>>::from(c)*tn))
            .reduce(|a, b| a.zip(b).map2(|(a, b)| a + b))
            .unwrap_or_default();

        let mut w = self*Self::FRAC_2_PI() + Self::one();
        let mut i = 0;
        while i < 4
        {
            w -= Self::one();
            if i % 2 == 0 && w < Self::zero()
            {
                w = -w;
            }
            w %= <Self as From<_>>::from(4.0);
            i += 1;
        }
        let two = <Self as From<_>>::from(2.0);
        let w = if w > Self::one() {two - w} else if w < -Self::one() {-two - w} else {w};

        let z = two*w*w - Self::one();

        p.polynomial(z)*w
    }

    fn tan(self) -> Self
    {
        let (sin, cos) = self.sin_cos();

        sin/cos
    }

    fn asin(self) -> Self
    {
        if self.is_nan()
        {
            return self
        }
        if self == Self::one()
        {
            return Self::FRAC_PI_2()
        }
        if self == -Self::one()
        {
            return -Self::FRAC_PI_2()
        }
        if self > Self::one() || self < -Self::one()
        {
            return Self::nan()
        }
        (self/(Self::one() - self*self).sqrt()).atan()
    }

    fn acos(self) -> Self
    {
        if self.is_nan()
        {
            return self
        }
        if self == Self::one()
        {
            return Self::zero()
        }
        if self == -Self::one()
        {
            return Self::PI()
        }
        if self > Self::one() || self < -Self::one()
        {
            return Self::nan()
        }
        ((Self::one() - self*self).sqrt()/self).atan()
    }

    fn atan(self) -> Self
    {
        if self.is_nan()
        {
            return self
        }
        if self.is_infinite()
        {
            return Self::FRAC_PI_2().copysign(self)
        }
        const TAYLOR: usize = 8;
        if self.abs() < Self::one()
        {
            let mut z = self;
            let mut y = z;

            for k in 1..TAYLOR
            {
                z *= -self*self;
                y += z/NumCast::from(1 + 2*k).unwrap()
            }

            y
        }
        else
        {
            let mut z = Self::one()/self;
            let mut y = Self::FRAC_PI_2() - z;

            for k in 1..TAYLOR
            {
                z /= -self*self;
                y -= z/NumCast::from(1 + 2*k).unwrap()
            }

            y
        }
    }

    fn atan2(self, other: Self) -> Self
    {
        if other.is_zero()
        {
            if self.is_zero()
            {
                return Self::zero()
            }
            if self.is_sign_negative()
            {
                return -Self::FRAC_PI_2()
            }
            else
            {
                return Self::FRAC_PI_2()
            }
        }
        let atan = (other/self).atan();
        if other.is_sign_positive()
        {
            return atan
        }
        return atan + if self.is_sign_positive()
        {
            Self::PI()
        }
        else
        {
            -Self::PI()
        }
    }

    fn sin_cos(self) -> (Self, Self)
    {
        const N: usize = 6;
        const C: [f64; N] = [
            1.276278962,
            -0.285261569,
            0.009118016,
            -0.000136587,
            0.000001185,
            -0.000000007
        ];
        let t: [[Self; N]; N] = ArrayOps::fill(
            |n| <[_; N]>::chebyshev_polynomial(1, n).unwrap()
        );
        let p: [Self; N] = t.zip(C)
            .map2(|(t, c)| t.map2(|tn| <Self as From<_>>::from(c)*tn))
            .reduce(|a, b| a.zip(b).map2(|(a, b)| a + b))
            .unwrap_or_default();

        let sin = {
            let mut w = self*Self::FRAC_2_PI();
            let mut i = 0;
            while i < 4
            {
                w -= Self::one();
                if i % 2 == 0 && w < Self::zero()
                {
                    w = -w;
                }
                w %= <Self as From<_>>::from(4.0);
                i += 1;
            }
            let two = <Self as From<_>>::from(2.0);
            let w = if w > Self::one() {two - w} else if w < -Self::one() {-two - w} else {w};
    
            let z = two*w*w - Self::one();
    
            p.polynomial(z)*w
        };

        let cos = {
            let mut w = self*Self::FRAC_2_PI() + Self::one();
            let mut i = 0;
            while i < 4
            {
                w -= Self::one();
                if i % 2 == 0 && w < Self::zero()
                {
                    w = -w;
                }
                w %= <Self as From<_>>::from(4.0);
                i += 1;
            }
            let two = <Self as From<_>>::from(2.0);
            let w = if w > Self::one() {two - w} else if w < -Self::one() {-two - w} else {w};

            let z = two*w*w - Self::one();

            p.polynomial(z)*w
        };

        (sin, cos)
    }

    fn exp_m1(self) -> Self
    {
        self.exp() - Self::one()
    }

    fn ln_1p(self) -> Self
    {
        (self + Self::one()).ln()
    }

    fn sinh(self) -> Self
    {
        if !self.is_finite()
        {
            return self
        }

        let ex = self.exp();
        let emx = (-self).exp();

        (ex - emx)*<Self as From<_>>::from(0.5)
    }

    fn cosh(self) -> Self
    {
        if !self.is_finite()
        {
            return self.abs()
        }

        let ex = self.exp();
        let emx = (-self).exp();

        (ex + emx)*<Self as From<_>>::from(0.5)
    }

    fn tanh(self) -> Self
    {
        if self.is_nan()
        {
            return self
        }

        let x2 = <Self as From<_>>::from(2.0)*self;
        let ex2 = x2.exp();
        let ex2p1 = ex2 + Self::one();
        let ex2m1 = ex2 - Self::one();

        if ex2m1.is_infinite()
        {
            return Self::one().copysign(self)
        }
        
        ex2m1/ex2p1
    }

    fn asinh(self) -> Self
    {
        if self.is_infinite() || self.is_nan()
        {
            return self
        }
        (self + (self*self + Self::one()).sqrt()).ln()
    }

    fn acosh(self) -> Self
    {
        if self.is_nan()
        {
            return self
        }
        if self < Self::one()
        {
            return Self::nan()
        }
        if self.is_infinite()
        {
            return self
        }
        (self + (self*self - Self::one()).sqrt()).ln()
    }

    fn atanh(self) -> Self
    {
        if self.is_nan()
        {
            return self
        }
        if self == Self::one()
        {
            return Self::infinity()
        }
        if self == -Self::one()
        {
            return Self::neg_infinity()
        }
        if self > Self::one() || self < -Self::one()
        {
            return Self::nan()
        }
        if self > Self::zero()
        {
            <Self as From<_>>::from(0.5)*((Self::one() + self)/(Self::one() - self)).ln()
        }
        else
        {
            -<Self as From<_>>::from(0.5)*((Self::one() - self)/(Self::one() + self)).ln()
        }
    }

    fn integer_decode(self) -> (u64, i16, i8)
    {
        let s = self.sign_bit();
        let e = self.exp_bits();
        let f = self.frac_bits();
        let bias = Self::exp_bias();

        let s = NumCast::from(s).unwrap();
        let e = NumCast::from(<i32 as NumCast>::from(e).unwrap() - <i32 as NumCast>::from(bias).unwrap()).unwrap();
        let f = NumCast::from(f).unwrap();

        (f, e, s)
    }
}

#[cfg(test)]
mod test
{
    use num_traits::{Float, Inv, One};

    use crate::fp::ieee754::{FpDouble, FpHalf, FpQuadruple, FpSingle};

    #[test]
    fn test_trig()
    {
        crate::tests::test_op1(Float::acos, Float::acos);
    }

    #[test]
    fn test_exp2()
    {
        crate::tests::test_op1(Float::exp2, Float::exp2);
    }

    #[test]
    fn test_log2()
    {
        crate::tests::test_op1(Float::log2, Float::log2);
    }

    #[test]
    fn test_powf()
    {
        crate::tests::test_op2(Float::powf, Float::powf);
    }

    #[test]
    fn test_round()
    {
        crate::tests::test_op1(Float::round, Float::round);
        crate::tests::test_op1(Float::ceil, Float::ceil);
        crate::tests::test_op1(Float::floor, Float::floor);
        crate::tests::test_op1(Float::trunc, Float::trunc);
    }

    #[test]
    fn test_sqrt()
    {
        let sqrt = FpDouble::from(3.0).sqrt();
        println!("{}^2 == {}", sqrt, sqrt*sqrt);
        let sqrt = (3.0).sqrt();
        println!("{}^2 == {}", sqrt, sqrt*sqrt);
    }
}