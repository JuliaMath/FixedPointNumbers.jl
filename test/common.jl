using FixedPointNumbers, Statistics, Random, Test
using FixedPointNumbers: bitwidth, rawtype, nbitsfrac
using Base.Checked

SP = VERSION >= v"1.6.0-DEV.771" ? " " : "" # JuliaLang/julia #37085

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
    f_max = bitwidth(T) - (T <: Signed)
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

# generator for cartesian product
function xypairs(::Type{X}) where X
    xs = typemin(X):eps(X):typemax(X)
    ((x, y) for x in xs, y in xs)
end


function test_floattype(TX::Type)
    @testset "floattype($X)" for X in target(TX, :i8, :i16, :i32, :i64; ex = :heavy)
        @test typemax(X) <= maxintfloat(floattype(X))
    end
end

function test_convert_from_nan(TX::Type)
    @testset "$X(nan)" for X in target(TX; ex = :thin)
        @test_throws ArgumentError X(Inf)
        @test_throws ArgumentError X(-Inf32)
        @test_throws ArgumentError X(NaN)
    end
end

function test_rem_type(TX::Type)
    @testset "% $X" for X in target(TX, :i8, :i16; ex = :thin)
        xs = typemin(X):0.1:typemax(X)
        @test all(x -> x % X === X(x), xs)
    end
end

function test_rem_nan(TX::Type)
    # TODO: avoid undefined behavior
    @testset "nan % $X" for X in target(TX, :i8, :i16, :i32, :i64; ex = :thin)
        @test NaN % X === NaN32 % X === NaN16 % X === zero(X)
    end
end

function test_neg(TX::Type)
    for X in target(TX, :i8; ex = :thin)
        xs = typemin(X):eps(X):typemax(X)
        fneg(x) = -float(x)
        @test all(x -> wrapping_neg(wrapping_neg(x)) === x, xs)
        @test all(x -> saturating_neg(x) === clamp(fneg(x), X), xs)
        @test all(x -> !(typemin(X) <= fneg(x) <= typemax(X)) ||
                       wrapping_neg(x) === checked_neg(x) === fneg(x) % X, xs)
    end
end

function test_abs(TX::Type)
    for X in target(TX, :i8; ex = :thin)
        xs = typemin(X):eps(X):typemax(X)
        fabs(x) = abs(float(x))
        @test all(x -> wrapping_abs(x) === (x > 0 ? x : wrapping_neg(x)), xs)
        @test all(x -> saturating_abs(x) === clamp(fabs(x), X), xs)
        @test all(x -> !(typemin(X) <= fabs(x) <= typemax(X)) ||
                       wrapping_abs(x) === checked_abs(x) === fabs(x) % X, xs)
    end
end

function test_add(TX::Type)
    for X in target(TX, :i8; ex = :thin)
        xys = xypairs(X)
        fadd(x, y) = float(x) + float(y)
        @test all(((x, y),) -> wrapping_sub(wrapping_add(x, y), y) === x, xys)
        @test all(((x, y),) -> saturating_add(x, y) === clamp(fadd(x, y), X), xys)
        @test all(((x, y),) -> !(typemin(X) <= fadd(x, y) <= typemax(X)) ||
                               wrapping_add(x, y) === checked_add(x, y) === fadd(x, y) % X, xys)
    end
end

function test_sub(TX::Type)
    for X in target(TX, :i8; ex = :thin)
        xys = xypairs(X)
        fsub(x, y) = float(x) - float(y)
        @test all(((x, y),) -> wrapping_add(wrapping_sub(x, y), y) === x, xys)
        @test all(((x, y),) -> saturating_sub(x, y) === clamp(fsub(x, y), X), xys)
        @test all(((x, y),) -> !(typemin(X) <= fsub(x, y) <= typemax(X)) ||
                               wrapping_sub(x, y) === checked_sub(x, y) === fsub(x, y) % X, xys)
    end
end

function test_mul(TX::Type)
    for X in target(TX, :i8; ex = :thin)
        xys = xypairs(X)
        fmul(x, y) = float(x) * float(y) # note that precision(Float32) < 32
        @test all(((x, y),) -> wrapping_mul(x, y) === fmul(x, y) % X, xys)
        @test all(((x, y),) -> saturating_mul(x, y) === clamp(fmul(x, y), X), xys)
        @test all(((x, y),) -> !(typemin(X) <= fmul(x, y) <= typemax(X)) ||
                               wrapping_mul(x, y) === checked_mul(x, y), xys)
    end
end

function test_isapprox(TX::Type)
    @testset "approx $X" for X in target(TX, :i8, :i16; ex = :light)
        xs = typemin(X):eps(X):typemax(X)-eps(X)
        @test all(x -> x ≈ x + eps(X), xs)
        @test all(x -> x + eps(X) ≈ x, xs)
        @test !any(x -> x - eps(X) ≈ x + eps(X), xs)
    end

end

function test_clamp_nan(TX::Type)
    @testset "clamp(nan, $X)" for X in target(TX; ex = :thin)
        @test clamp( Inf, X) === clamp( Inf32, X) === typemax(X)
        @test clamp(-Inf, X) === clamp(-Inf32, X) === typemin(X)
        @test clamp( NaN, X) === clamp( NaN32, X) === zero(X)
    end
end

function test_isinteger(TX::Type)
    @testset "isinteger(::$X)" for X in target(TX, :i8, :i16)
        xs = typemin(X):eps(X):typemax(X)
        @test all(x -> isinteger(x) == isinteger(float(x)), xs)
    end
end

function test_rand(TX::Type)
    @testset "rand(::$X)" for X in target(TX; ex = :thin)
        @test isa(rand(X), X)
        a = rand(X, (3, 5))
        @test ndims(a) == 2 && eltype(a) === X
        @test size(a) == (3, 5)
    end
end
