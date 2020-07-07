using FixedPointNumbers, Statistics, Random, Test
using FixedPointNumbers: bitwidth, rawtype, nbitsfrac
using Base.Checked

"""
    target(X::Type, Ss...; ex = :default)

Return a generator which enumerates the target types for testing.

# Arguments
- `X`: target base type
- `Ss`: symbols for specifying the target raw types
  - `:i*` : a `Signed` type if `X === Fixed`, or an `Unsigned` type if `X === Normed`
  - `:s*` : a `Signed` type (not yet supported)
  - `:u*` : an `Unsigned` type (not yet supported)
- `ex`: exhaustivity of `f`s (see also the [`target_f`](@ref) function)
  - `:heavy`: all supported `f`s
  - `:default`: same as `:heavy` for 8-/16-bit types, and same as `:light` otherwise
  - `:light`: important `f`s for byte boundaries and floating point types
  - `:thin`: maximum and half `f`s per type

# Example
```julia
julia> collect(target(Normed, :i8, :i32; ex = :default))
21-element Array{DataType,1}:
 Normed{UInt8,1}
 Normed{UInt8,2}
 Normed{UInt8,3}
 Normed{UInt8,4}
 Normed{UInt8,5}
 Normed{UInt8,6}
 Normed{UInt8,7}
 Normed{UInt8,8}
 Normed{UInt32,1}
 Normed{UInt32,7}
 Normed{UInt32,8}
 Normed{UInt32,9}
 Normed{UInt32,10}
 Normed{UInt32,11}
 Normed{UInt32,15}
 Normed{UInt32,16}
 Normed{UInt32,17}
 Normed{UInt32,23}
 Normed{UInt32,24}
 Normed{UInt32,31}
 Normed{UInt32,32}
```
"""
function target(X::Type, Ss...; ex = :default)
    Ts = symbol_to_inttype.(X, Ss)
    (X{T,f} for T in Ts for f in target_f(X, T; ex = ex))
end
target(X::Type; ex = :default) = target(X, :i8, :i16, :i32, :i64, :i128; ex = ex)

"""
    target_f(X::Type, T::Type; ex = :default)

Return a tuple or range of the number of fractional bits `f` to be tested.

# Arguments
The `X` specifies the target base type, i.e. `Fixed` or `Normed`, and the `T`
specifies the target raw type.

## `ex` keyword
The `ex` specifies the exhaustivity of `f`s.
The following are examples of `target_f(Normed, T)`. The marker `x` means the
target and the marker `-` means not the target.

### `:heavy` -- all supported `f`s
```
            |    3                   2                   1                  |
          f |2 1 0 9 8 7 6 5:4 3 2 1 0 9 8 7:6 5 4 3 2 1 0 9:8 7 6 5 4 3 2 1|
T == UInt8  |               :               :               :x x x x x x x x|
T == UInt16 |               :               :x x x x x x x x:x x x x x x x x|
T == UInt32 |x x x x x x x x:x x x x x x x x:x x x x x x x x:x x x x x x x x|
```
### `:default` -- same as `:heavy` for 8-/16-bit types, and same as `:light` otherwise
```
            |    3                   2                   1                  |
          f |2 1 0 9 8 7 6 5:4 3 2 1 0 9 8 7:6 5 4 3 2 1 0 9:8 7 6 5 4 3 2 1|
T == UInt8  |               :               :               :x x x x x x x x|
T == UInt16 |               :               :x x x x x x x x:x x x x x x x x|
T == UInt32 |x x - - - - - -:x x - - - - - x:x x - - - x x x:x x - - - - - x|
```

### `:light` -- important `f`s for byte boundaries and floating point types
```
            |    3                   2                   1                  |
          f |2 1 0 9 8 7 6 5:4 3 2 1 0 9 8 7:6 5 4 3 2 1 0 9:8 7 6 5 4 3 2 1|
T == UInt8  |               :               :               :x x - - - - - x|
T == UInt16 |               :               :x x - - - x x x:x x - - - - - x|
T == UInt32 |x x - - - - - -:x x - - - - - x:x x - - - x x x:x x - - - - - x|
                             |                         |
                             +--precision(Float32)     +--precision(Float16)
```

### `:thin` -- maximum and half `f`s per type
```
            |    3                   2                   1                  |
          f |2 1 0 9 8 7 6 5:4 3 2 1 0 9 8 7:6 5 4 3 2 1 0 9:8 7 6 5 4 3 2 1|
T == UInt8  |               :               :               :x - - - x - - -|
T == UInt16 |               :               :x - - - - - - -:x - - - - - - -|
T == UInt32 |x - - - - - - -:- - - - - - - -:x - - - - - - -:- - - - - - - -|
```
"""
function target_f(X::Type, T::Type{<:Integer}; ex = :default)
    f_min = X === Fixed ? 0 : 1
    f_max = bitwidth(T) - (T <: Signed) - 1 + f_min
    ex === :heavy && return f_min:f_max
    ex === :default && bitwidth(T) <= 16 && return f_min:f_max
    ex === :thin && return ((f_max + 1) ÷ 2, f_max)
    if ex === :light || ex === :default
        itr = Iterators.filter(x -> x <= f_max, target_f_series(X, T))
        return (itr...,)
    end
    error()
end

target_f_series(::Type{Fixed}, T::Type{<:Integer}) =
    (0, 1, 7, 8, 9,
     10, 11, 15, 16, 17,
     23, 24, 31, 32, 33,
     52, 53, 63, 64, 65,
     112, 113, 127)

target_f_series(::Type{Normed}, T::Type{<:Integer}) =
    (1, 7, 8, 9,
     10, 11, 15, 16, 17,
     23, 24, 31, 32, 33,
     52, 53, 63, 64, 65,
     112, 113, 127, 128)
