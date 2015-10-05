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

# Type hierarchy
This library defines an abstract type `FixedPoint{T <: Integer, f}` as a subtype of `Real`. The parameter `T` is the underlying representation and `f` is the number of fraction bits.

For signed integers, there is a fixed-point type `Fixed{T, f}` and for unsigned integers, there is the `UFixed{T, f}` type.

These types, built with `f` fraction bits, map the closed interval [0.0,1.0]
to the span of numbers with `f` bits.
For example, the `UFixed8` type is represented internally by a `UInt8`, and makes
`0x00` equivalent to `0.0` and `0xff` to `1.0`.
The types `UFixed10`, `UFixed12`, `UFixed14`, and `UFixed16` are all based on `UInt16`
and reach the value `1.0` at 10, 12, 14, and 16 bits, respectively (`0x03ff`, `0x0fff`,
`0x3fff`, and `0xffff`).

To construct such a number, use `convert(UFixed12, 1.3)`, `ufixed12(1.3)`, or the literal syntax `0x14ccuf12`.
The latter syntax means to construct a `UFixed12` (it ends in `uf12`) from the `UInt16` value
`0x14cc`.

There currently is no literal syntax for signed `Fixed` numbers. 

[wikipedia]: http://en.wikipedia.org/wiki/Fixed-point_arithmetic
