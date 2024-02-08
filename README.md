# Custom Float

This crate adds a custom floating point number type, Fp<U, EXP_SIZE, FRAC_SIZE>, where the bit size of the exponent and mantissa can be set separately.

The number then automatically implements num::Float, and supports all ordinary floating point operations.

This allows simple implementation of special floating point types, such as TensorFloat, IEEE754 quadruple precision float, and BFloat16.

The accuracy of the various floating point operation is not perfect, but works well enough to be usable.

All arithmetic operations are performed using bitwise operators and integer arithmetic.
