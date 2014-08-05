# FixedPointNumbers

This library exports fixed-point number types.
A [fixed-point number][wikipedia] represents a fractional, or non-integral, number.
In contrast with the more widely known floating-point numbers, fixed-point
numbers have a fixed number of digits (bits) after the decimal (radix) point.
They are effectively integers scaled by a constant factor.

Fixed-point numbers can be used to perform arithmetic. Another practical
application is to implicitly rescale integers without modifying the
underlying representation.

This library exports two categories of fixed-point types. Fixed-point types are
used like any other number: they can be added, multiplied, raised to a power,
etc. In many cases these operations result in conversion to floating-point types.

## Fixed32 (signed fixed-point numbers)

For signed integers, there is a 32-bit fixed-point type `Fixed32{f}`.
The parameter `f` is the number of fraction bits. There is also an abstract subtype of
`Real` called `Fixed`.

To use it, convert numbers to a `Fixed32` type, or call `Fixed32(x)`, which will default
to constructing a `Fixed32{16}`.

## Ufixed (unsigned fixed-point numbers)

For unsigned integers, there is a family of subtypes of the abstract `Ufixed` type.
These types, built with `f` fraction bits, map the closed interval [0.0,1.0]
to the span of numbers with `f` bits.
For example, the `Ufixed8` type is represented internally by a `Uint8`, and makes
`0x00` equivalent to `0.0` and `0xff` to `1.0`.
The types `Ufixed10`, `Ufixed12`, `Ufixed14`, and `Ufixed16` are all based on `Uint16`
and reach the value `1.0` at 10, 12, 14, and 16 bits, respectively (`0x03ff`, `0x0fff`,
`0x3fff`, and `0xffff`).

To construct such a number, use `convert(Ufixed12, 1.3)` or the literal syntax `0x14ccuf12`.
The latter syntax means to construct a `Ufixed12` (it ends in `uf12`) from the `Uint16` value
`0x14cc`.

[wikipedia]: http://en.wikipedia.org/wiki/Fixed-point_arithmetic
