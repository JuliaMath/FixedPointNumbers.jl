# FixedPointNumbers

[![Build Status](https://travis-ci.org/JuliaMath/FixedPointNumbers.jl.svg?branch=master)](https://travis-ci.org/JuliaMath/FixedPointNumbers.jl)

[![codecov.io](http://codecov.io/github/JuliaMath/FixedPointNumbers.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaMath/FixedPointNumbers.jl?branch=master)

This library implements fixed-point number types.  A
[fixed-point number][wikipedia] represents a fractional, or
non-integral, number.  In contrast with the more widely known
floating-point numbers, with fixed-point numbers the decimal point
doesn't "float": fixed-point numbers are effectively integers that are
interpreted as being scaled by a constant factor.  Consequently, they
have a fixed number of digits (bits) after the decimal (radix) point.

Fixed-point numbers can be used to perform arithmetic. Another practical
application is to implicitly rescale integers without modifying the
underlying representation.

This library exports two categories of fixed-point types. Fixed-point types are
used like any other number: they can be added, multiplied, raised to a power,
etc. In some cases these operations result in conversion to floating-point types.

# Type hierarchy and interpretation

This library defines an abstract type `FixedPoint{T <: Integer, f}` as a
subtype of `Real`. The parameter `T` is the underlying machine representation and `f`
is the number of fraction bits.

For `T<:Signed` (a signed integer), there is a fixed-point type
`Fixed{T, f}`; for `T<:Unsigned` (an unsigned integer), there is the
`Normed{T, f}` type. However, there are slight differences in behavior
that go beyond signed/unsigned distinctions.

The `Fixed{T,f}` types use 1 bit for sign, and `f` bits to represent
the fraction. For example, `Fixed{Int8,7}` uses 7 bits (all bits
except the sign bit) for the fractional part. The value of the number
is interpreted as if the integer representation has been divided by
`2^f`. Consequently, `Fixed{Int8,7}` numbers `x` satisfy

```
-1.0 = -128/128 ≤ x ≤ 127/128 ≈ 0.992.
```

because the range of `Int8` is from -128 to 127.

In contrast, the `Normed{T,f}`, with `f` fraction bits, map the closed
interval [0.0,1.0] to the span of numbers with `f` bits.  For example,
the `N0f8` type (aliased to `Normed{UInt8,8}`) is represented
internally by a `UInt8`, and makes `0x00` equivalent to `0.0` and
`0xff` to `1.0`. Consequently, `Normed` numbers are scaled by `2^f-1`
rather than `2^f`.  The type aliases `N6f10`, `N4f12`,
`N2f14`, and `N0f16` are all based on `UInt16` and reach the
value `1.0` at 10, 12, 14, and 16 bits, respectively (`0x03ff`,
`0x0fff`, `0x3fff`, and `0xffff`). The `NXfY` notation is used for
compact printing and the `fY` component informs about the number of
fractional bits and `X+Y` equals the number of underlying bits used.

To construct such a number, use `convert(N4f12, 1.3)`, `N4f12(1.3)`,
`Normed{UInt16,12}(1.3)`, or `reinterpret(N4f12, 0x14cc)`.
The latter syntax means to construct a `N4f12` (it ends in
`uf12`) from the `UInt16` value `0x14cc`.

More generally, an arbitrary number of bits from any of the standard unsigned
integer widths can be used for the fractional part.  For example:
`Normed{UInt32,16}`, `Normed{UInt64,3}`, `Normed{UInt128,7}`.

[wikipedia]: http://en.wikipedia.org/wiki/Fixed-point_arithmetic
