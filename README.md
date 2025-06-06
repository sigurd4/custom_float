[![Build Status (nightly)](https://github.com/sigurd4/custom_float/workflows/Build-nightly/badge.svg)](https://github.com/sigurd4/custom_float/actions/workflows/build-nightly.yml)
[![Build Status (nightly, all features)](https://github.com/sigurd4/custom_float/workflows/Build-nightly-all-features/badge.svg)](https://github.com/sigurd4/custom_float/actions/workflows/build-nightly-all-features.yml)

[![Build Status (stable)](https://github.com/sigurd4/custom_float/workflows/Build-stable/badge.svg)](https://github.com/sigurd4/custom_float/actions/workflows/build-stable.yml)
[![Build Status (stable, all features)](https://github.com/sigurd4/custom_float/workflows/Build-stable-all-features/badge.svg)](https://github.com/sigurd4/custom_float/actions/workflows/build-stable-all-features.yml)

[![Test Status](https://github.com/sigurd4/custom_float/workflows/Test/badge.svg)](https://github.com/sigurd4/custom_float/actions/workflows/test.yml)
[![Lint Status](https://github.com/sigurd4/custom_float/workflows/Lint/badge.svg)](https://github.com/sigurd4/custom_float/actions/workflows/lint.yml)

[![Latest Version](https://img.shields.io/crates/v/custom_float.svg)](https://crates.io/crates/custom_float)
[![License:MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Documentation](https://img.shields.io/docsrs/custom_float)](https://docs.rs/custom_float)
[![Coverage Status](https://img.shields.io/codecov/c/github/sigurd4/custom_float)](https://app.codecov.io/github/sigurd4/custom_float)

# custom_float

This crate adds a custom floating point number type, `Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>`, where the bit size of the exponent and mantissa can be set separately, as well as the base of the exponent (which is normally 2).

This allows simple implementation of special floating point types, such as TensorFloat, IEEE754 Decimals (decimal32, decimal64, etc.), Fp80, and BFloat16.

## Composition

`U` is the underlying unsigned integer type which is used to represent the number.

`SIGN_BIT` is wether or not the number has a sign bit.

`EXP_SIZE` is the size of the exponent in bits.

`INT_SIZE` is the size of the integer part of the mantissa in bits. If zero, then the integer bit is implicit.

`FRAC_SIZE` is the size of the fractional part of the mantissa in bits.

`EXP_BASE` is the base of the exponent.

The total bit size of `U` must be greater or equal to `SIGN_BIT` + `EXP_SIZE` + `INT_SIZE` + `FRAC_SIZE` to contain the entire number.

The bit layout is as follows:
```txt
No data: | Sign:      | Exponent:  | Integer:   | Fractional: |
<  ..  > | <SIGN_BIT> | <EXP_SIZE> | <INT_SIZE> | <FRAC_SIZE> |
```

The value of a real floating-point number is the following:
```txt
x = (-1)**sign*EXP_BASE**(exponent - bias)*mantissa
```

where the bias equals
```txt
bias = 2**(EXP_SIZE - 1) - 1
```

If the exponent has the maximum value, the number is either infinity or NaN.

## Features

This crate provides the type `Fp`, and not really anything else.

### Traits

All `Fp`'s automatically implement `num::Float`, and supports all ordinary floating point operations you'd expect (`sin`, `cos`, `tanh`, `sqrt`, `exp`, `ln`, `powf`, etc. as well as operators `+`, `-`, `*`, `/`, `%`). I've also implemented some equivalents to some of the rational functions from `libm` (like `erf`, `erfc`, and bessel-functions `j0`, `y0`, etc.).

### Size

The biggest integers in the standard library are 128-bit, so if you want floating point numbers with more bits than that, you have to provide your own 128-bit unsigned integer type, but in theory, it should be possible to have, say, 16384-bit floats, if you want to.

### Conversion

All `Fp`'s can be converted into each other with `Fp::from_fp`. For now, due to trait-implementation conflicts, the `From` trait can't be used for this, because it would conflict with the implementation of `T: From<T>`.

All `Fp`'s implement `From` and `Into` for all standard-library numeric types (unsigned integers: `u8`, `u16`, `u32`, `u64`, `u128`, signed integers `i8`, `i16`, `i32`, `i64`, `i128`, and floats `f16`, `f32`, `f64`, `f128`).

Of course, narrowing conversions will result in rounding-errors and unbounded values.

### Aliases

There are a couple of pre-made aliases for floats available too. They're just type-aliases for `Fp`.

### Compatability

The following types are bitwise interchangeable:
- `f16` and `FpHalf`
- `f32` and `FpSingle`
- `f64` and `FpDouble`

I've read that `binary128` may have an odd two-bit prefix in the exponent, that may work slightly different from the lesser IEEE754 binary floats. According to my tests, `f128` and `FpQuadruple` also seem to be interchangeable, even though i'm not sure if `f128`s always work like this on every compilation target. Please make a report if you notice problems with this!

## Examples

```rust
use custom_float::Fp;

// This type is also available as `custom_float::ieee754::FpSingle`
type FpSingle = Fp<u32, true, 8, 0, 23, 2>;

let two = FpSingle::from(2);
let four = FpSingle::from(4);

assert_eq!(two + two, four);
assert_eq!(two * two, four);
assert_eq!(two / four, two.recip());
assert_eq!(four - two, two);
```

## Performance

This crate is obviously slower than the primitive floats in the standard library, because my implementations are not just LLVM-intrinisics. For nonstandard floating-point formats, you probably wont even find processors that can do these kinds of operations atomically. For base 2 floats that also match the primitive standard-library floats (f16, f32, f64, f128) in form, you can enable the feature `use_std_float` to convert them and do the operation natively instead, which may give a performance boost. This will only happen with floats that are trivially convertable with either f16, f32, f64 or f128.

You can see comparisons with this library's floats and the standard library's in the [plots/bench](https://github.com/sigurd4/custom_float/tree/master/plots/bench) folder. For some methods not found in `std`/`core`, i use `libm` instead as a comparison.

## Accuracy

The accuracy is obviously dependent on the size of your float's mantissa and its base, and its range depends on the size of its exponent. But given a float with a (practically) infinite resolution, some functions on here still give an error. My goal is to make that error as small as possible. The accuracy of all of the floating point operations are not perfect, but work well enough to be usable.

You can see comparisons with this library's floats and the standard library's in the [plots](https://github.com/sigurd4/custom_float/tree/master/plots) folder, and error in the [plots/error](https://github.com/sigurd4/custom_float/tree/master/plots/error). For some methods not found in `std`/`core`, i use `libm` instead as a comparison.

## Planned features

- Make more and more of the functions work at compile-time.
    - This is currently difficult because my code is very generic, and relies on traits from the `num-traits` crate (that are not const-traits, because those are experimental). Once rust's const-traits are a stable language-feature, you'll see more of this.
- Stabilize large bases. Currently i mostly just test base 2 and 10. For large bases (say, 1000) you tend to get integer overflow.
- Stabilize edge cases like unsigned floats, exponentless floats (fixed-point, really), and mantissaless floats (just exponentials?).
- Bigfloats? (maybe for a separate crate)
- Serde-integration
- Proper (not convert-to-nearest-std-float) implementations of `Debug`, `Display`, and `FromStr`.
- Use signaling NaN's correctly.
- Make it run faster, of course!

Suggestions are welcome.