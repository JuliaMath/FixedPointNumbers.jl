# FixedPoint

This library exports a 32-bit fixed-point type `Fixed32{f}`.
The parameter `f` is the number of fraction bits. There is also an abstract subtype of
`Real` called `Fixed`.

To use it, convert numbers to a `Fixed32` type, or call `Fixed32(x)`, which will default
to constructing a `Fixed32{16}`. Then use ordinary arithmetic.
