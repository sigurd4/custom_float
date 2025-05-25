# custom_float

This crate adds a custom floating point number type, `Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>`, where the bit size of the exponent and mantissa can be set separately.

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

The number then automatically implements `num::Float`, and supports all ordinary floating point operations.

This allows simple implementation of special floating point types, such as TensorFloat, IEEE754 Quadruple/binary128, Fp80, and BFloat16.

The accuracy of all of the floating point operations are not perfect, but work well enough to be usable. Various plots showing the accuracy of basic functions are shown in the [plots](https://github.com/sigurd4/custom_float/tree/master/plots) subfolder.

All floats can be converted into each other painlessly, though the conversion may produce rounding errors or unbounded outputs when converting to a float with lesser resolution.

# Examples

```rust
#![feature(generic_const_exprs)]

use custom_float::Fp;

type FpSingle = Fp<u32, true, 8, 0, 23, 2>;

let two = FpSingle::from(2);
let four = FpSingle::from(4);

assert_eq!(two + two, four);
```