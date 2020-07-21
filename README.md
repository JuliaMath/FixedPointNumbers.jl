# FixedPointNumbers

[![Build Status][action-img]][action-url]
[![Build Status][pkgeval-img]][pkgeval-url]
[![coverage][codecov-img]][codecov-url]

This library implements fixed-point number types.  A
[fixed-point number] represents a fractional, or
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

```julia
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

To construct such a number, use `1.3N4f12`, `N4f12(1.3)`, `convert(N4f12, 1.3)`,
`Normed{UInt16,12}(1.3)`, or `reinterpret(N4f12, 0x14cc)`.
The last syntax means to construct a `N4f12` from the `UInt16` value `0x14cc`.

More generally, an arbitrary number of bits from any of the standard unsigned
integer widths can be used for the fractional part.  For example:
`Normed{UInt32,16}`, `Normed{UInt64,3}`, `Normed{UInt128,7}`.

# Computation with Fixed and Normed numbers

You can perform mathematical operations with `FixedPoint` numbers, but keep in mind
that they are vulnerable to both [rounding] and [overflow]. For example:

```julia
julia> x = N0f8(0.8)
0.8N0f8

julia> float(x) + x
1.6f0

julia> x + x
0.596N0f8
```

This is a consequence of the rules that govern overflow in integer arithmetic:

```julia
julia> y = reinterpret(x)        # `reinterpret(x::FixedPoint)` reinterprets as the underlying "raw" type
0xcc

julia> reinterpret(N0f8, y + y)  # add two UInt8s and then reinterpret as N0f8
0.596N0f8
```

Similarly,

```julia
julia> x = eps(N0f8)             # smallest nonzero `N0f8` number
0.004N0f8

julia> x*x
0.0N0f8
```

which is rounding-induced [underflow].  Finally,

```julia
julia> x = N4f12(15)
15.0N4f12

julia> x*x
ERROR: ArgumentError: Normed{UInt16,12} is a 16-bit type representing 65536 values from 0.0 to 16.0037; cannot represent 225.0
Stacktrace:
 [1] throw_converterror(::Type{Normed{UInt16,12}}, ::Float32) at /home/tim/.julia/dev/FixedPointNumbers/src/FixedPointNumbers.jl:251
 [2] _convert at /home/tim/.julia/dev/FixedPointNumbers/src/normed.jl:77 [inlined]
 [3] FixedPoint at /home/tim/.julia/dev/FixedPointNumbers/src/FixedPointNumbers.jl:51 [inlined]
 [4] convert at ./number.jl:7 [inlined]
 [5] *(::Normed{UInt16,12}, ::Normed{UInt16,12}) at /home/tim/.julia/dev/FixedPointNumbers/src/normed.jl:254
 [6] top-level scope at REPL[16]:1
```

In some circumstances, it may make most sense to think of `FixedPoint` numbers as *storage types*
rather than computational types. You can call `float(x)` to convert `x` to a floating-point equivalent that is reasonably
safe for computation; in the type domain, `floattype(T::Type)` returns the corresponding type.
Note that in some cases `floattype(T)` differs from `float`'s behavior on the corresponding "raw" type:

```julia
julia> float(UInt8)
Float64

julia> floattype(N0f8)
Float32
```

Because of the role of FixedPointNumbers in domains such as image-processing, this package tries to limit the expansion of the
number of bits needed to store results.


## Contributing to this package

Please see [CONTRIBUTING.md](CONTRIBUTING.md) for information about improving this package.


[fixed-point number]: http://en.wikipedia.org/wiki/Fixed-point_arithmetic
[overflow]: https://en.wikipedia.org/wiki/Integer_overflow
[rounding]: https://en.wikipedia.org/wiki/Round-off_error
[underflow]: https://en.wikipedia.org/wiki/Arithmetic_underflow


<!-- badges -->

[action-img]: https://github.com/JuliaMath/FixedPointNumbers.jl/workflows/Unit%20test/badge.svg
[action-url]: https://github.com/JuliaMath/FixedPointNumbers.jl/actions

[pkgeval-img]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/F/FixedPointNumbers.svg
[pkgeval-url]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html

[codecov-img]: https://codecov.io/gh/JuliaMath/FixedPointNumbers.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/JuliaMath/FixedPointNumbers.jl

