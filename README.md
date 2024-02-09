# Custom Float

This crate adds a custom floating point number type, `Fp<U, EXP_SIZE, INT_BIT, FRAC_SIZE>`, where the bit size of the exponent and mantissa can be set separately.

`EXP_SIZE` is the size of the exponent in bits.

`INT_BIT` is wether or not the int bit is explicit (`true`), or implicit (`false`).

`FRAC_SIZE` is the size of the fractional part/mantissa in bits.

The number then automatically implements `num::Float`, and supports all ordinary floating point operations.

This allows simple implementation of special floating point types, such as TensorFloat, IEEE754 Quadruple/binary128, Fp80, and BFloat16.

The accuracy of all of the floating point operations are not perfect, but works well enough to be usable.

All arithmetic operations are performed using bitwise operators and integer arithmetic.

All floats can be converted into each other efficiently.
