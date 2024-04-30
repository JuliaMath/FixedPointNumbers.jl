module FixedPointNumbers

import Base: ==, <, <=, -, +, *, /, ~, isapprox,
             convert, promote_rule, print, show, bitstring, abs, decompose,
             isnan, isinf, isfinite, isinteger,
             zero, oneunit, one, typemin, typemax, floatmin, floatmax, eps, reinterpret,
             big, rationalize, float, trunc, round, floor, ceil, bswap, clamp,
             div, fld, cld, rem, mod, mod1, fld1, min, max, minmax,
             signed, unsigned, copysign, flipsign, signbit,
             length

import Random: Random, AbstractRNG, SamplerType, rand!

using Base: @pure

"""
    FixedPoint{T <: Integer, f} <: Real

Supertype of the two fixed-point number types: `Fixed{T, f}` and `Normed{T, f}`.

The parameter `T` is the underlying machine representation and `f` is the number
of fraction bits.
"""
abstract type FixedPoint{T <: Integer, f} <: Real end


export
    FixedPoint,
    Fixed,
    Normed,
    floattype,
# "special" typealiases
    # Q and N typealiases are exported in separate source files
# Functions
    scaledual

include("utilities.jl")
using .Utilities
import .Utilities: floattype, rawone, nbitsfrac, rawtype, signbits, nbitsint, scaledual

# reinterpretation
reinterpret(x::FixedPoint) = x.i
reinterpret(::Type{T}, x::FixedPoint{T,f}) where {T,f} = x.i
reinterpret(::Type{X}, x::T) where {T <: Integer, X <: FixedPoint{T}} = X(x, 0)

# static parameters
nbitsfrac(::Type{X}) where {T, f, X <: FixedPoint{T,f}} = f
rawtype(::Type{X}) where {T, X <: FixedPoint{T}} = T

# traits based on static parameters
signbits(::Type{X}) where {T, X <: FixedPoint{T}} = T <: Unsigned ? 0 : 1
nbitsint(::Type{X}) where {X <: FixedPoint} = bitwidth(X) - nbitsfrac(X) - signbits(X)

# constructor-style conversions
(::Type{X})(x::X) where {X <: FixedPoint}      = x
(::Type{X})(x::Number) where {X <: FixedPoint} = _convert(X, x)

function (::Type{<:FixedPoint})(x::AbstractChar)
    throw(ArgumentError("FixedPoint (Fixed or Normed) cannot be constructed from a Char"))
end
(::Type{X})(x::Complex) where {X <: FixedPoint} = X(convert(real(typeof(x)), x))
function (::Type{X})(x::Base.TwicePrecision) where {X <: FixedPoint}
    floattype(X) === BigFloat ? X(big(x)) : X(convert(floattype(X), x))
end

# conversions
function Base.Bool(x::FixedPoint)
    x == zero(x) ? false : x == oneunit(x) ? true : throw(InexactError(:Bool, Bool, x))
end
function (::Type{Ti})(x::FixedPoint) where {Ti <: Integer}
    isinteger(x) || throw(InexactError(:Integer, typeof(x), x))
    floor(Ti, x)
end
Base.Rational{Ti}(x::FixedPoint) where {Ti <: Integer} = Rational{Ti}(Rational(x))

big(::Type{<:FixedPoint}) = BigFloat
big(x::FixedPoint) = convert(BigFloat, x)

rationalize(x::FixedPoint; tol::Real=eps(x)) = rationalize(Int, x, tol=tol)
function rationalize(::Type{Ti}, x::FixedPoint; tol::Real=eps(x)) where Ti <: Integer
    tol <= eps(x) ? Rational{Ti}(x) : rationalize(Ti, float(x), tol)
end

"""
    isapprox(x::FixedPoint, y::FixedPoint; rtol=0, atol=max(eps(x), eps(y)))

For FixedPoint numbers, the default criterion is that `x` and `y` differ by no
more than `eps`, the separation between adjacent fixed-point numbers.
"""
@inline function isapprox(x::X, y::X; rtol=0, atol=eps(X)) where {X <: FixedPoint}
    n, m = minmax(x, y) # m >= n
    if rtol == zero(rtol)
        _isapprox_atol(m, n, atol)
    elseif atol == zero(atol)
        _isapprox_rtol(m, n, rtol)
    else
        _isapprox_atol(m, n, atol) | _isapprox_rtol(m, n, rtol)
    end
end
function isapprox(x::FixedPoint, y::FixedPoint; rtol=0, atol=max(eps(x), eps(y)))
    isapprox(promote(x, y)...; rtol=rtol, atol=atol)
end
function _isapprox_atol(m::X, n::X, atol::X) where {X <: FixedPoint}
    unsigned(m.i - n.i) <= unsigned(max(atol, zero(X)).i)
end
function _isapprox_atol(m::X, n::X, atol) where {X <: FixedPoint}
    unsigned(m.i - n.i) <= div(atol, eps(X))
end
function _isapprox_rtol(m::X, n::X, rtol) where {X <: FixedPoint}
    unsigned(m.i - n.i) <= rtol * max(unsigned(abs(m.i)), unsigned(abs(n.i)))
end

# predicates
isinteger(x::FixedPoint) = x == trunc(x) # TODO: use floor(x) when dropping support for Fixed{Int8,8}
isfinite(x::FixedPoint) = true
isnan(x::FixedPoint) = false
isinf(x::FixedPoint) = false

# identities
zero(::Type{X}) where {X <: FixedPoint} = X(zero(rawtype(X)), 0)
oneunit(::Type{X}) where {X <: FixedPoint} = X(rawone(X), 0)
one(::Type{X}) where {X <: FixedPoint} = oneunit(X)

# traits
eps(::Type{X}) where {X <: FixedPoint} = X(oneunit(rawtype(X)), 0)
typemax(::Type{T}) where {T <: FixedPoint} = T(typemax(rawtype(T)), 0)
typemin(::Type{T}) where {T <: FixedPoint} = T(typemin(rawtype(T)), 0)
floatmin(::Type{T}) where {T <: FixedPoint} = eps(T)
floatmax(::Type{T}) where {T <: FixedPoint} = typemax(T)


"""
    floattype(::Type{T})::Type{<:AbstractFloat}

Return a minimal type suitable for performing computations with instances of type `T` without integer overflow.

The fallback definition of `floattype(T)` applies only to `T<:AbstractFloat`.
However, it is permissible to extend `floattype` to return types that are not subtypes of
`AbstractFloat`; the key characteristic is that the return type should support computation without integer overflow.

In general the returned type should have the minimum bitwidth needed to encode the full precision of the input type.
however, a priority should be placed on computational efficiency; consequently, types like `Float16` should be avoided
except in scenarios where they are guaranteed to have hardware support.

# Examples

A classic usage is to avoid overflow behavior by promoting `FixedPoint` to `AbstractFloat`

```jldoctest; setup = :(using FixedPointNumbers)
julia> x = N0f8(1.0)
1.0N0f8

julia> x + x # overflow
0.996N0f8

julia> T = floattype(x)
Float32

julia> T(x) + T(x)
2.0f0
```

The following represents a valid extension of `floattype` to non-AbstractFloats:

```julia
julia> using FixedPointNumbers, ColorTypes

julia> floattype(RGB{N0f8})
RGB{Float32}
```

`RGB` itself is not a subtype of `AbstractFloat`, but unlike `RGB{N0f8}` operations with `RGB{Float32}` are not subject to integer overflow.
"""
floattype(::Type{X}) where {T <: ShortInts, X <: FixedPoint{T}} = Float32
floattype(::Type{X}) where {T <: Integer, X <: FixedPoint{T}} = Float64
floattype(::Type{X}) where {T <: LongInts, X <: FixedPoint{T}} = BigFloat

float(x::FixedPoint) = convert(floattype(x), x)

function minmax(x::X, y::X) where {X <: FixedPoint}
    a, b = minmax(reinterpret(x), reinterpret(y))
    X(a,0), X(b,0)
end

bitstring(x::FixedPoint) = bitstring(x.i)

bswap(x::X) where {X <: FixedPoint} = sizeof(X) == 1 ? x : X(bswap(x.i), 0)

# At least on Julia v1.5.0 or earlier, the following specialization helps the
# SIMD vectorization. (cf. PR #194)
clamp(x::X, lo::X, hi::X) where {X <: FixedPoint} = X(clamp(x.i, lo.i, hi.i), 0)

clamp(x, ::Type{X}) where {X <: FixedPoint} = clamp(x, typemin(X), typemax(X)) % X

# Since `FixedPoint` is not an integer type, it is not clear in what type
# `signed` and `unsigned` for `FixedPoint` should return values. They should
# currently throw errors in case we support "unsigned Fixed" or "signed Normed"
# in the future. The following "incomplete" code is necessary for Julia v1.0
# etc. to prevent accidental conversion to an integer type.
signed(x::X) where {X <: FixedPoint} = signed(X)(signed(x.i), 0)
unsigned(x::X) where {X <: FixedPoint} = unsigned(X)(unsigned(x.i), 0)

function copysign(x::X, y::Real) where {T, X <: FixedPoint{T}}
    T <: Signed ? X(copysign(x.i, y), 0) : throw_not_a_signed_number_error(x)
end
function flipsign(x::X, y::Real) where {T, X <: FixedPoint{T}}
    T <: Signed ? X(flipsign(x.i, y), 0) : throw_not_a_signed_number_error(x)
end
if copysign(-1, 0x1) !== 1 # for Julia v1.0 and v1.1 (julia #30748)
    copysign(x::X, y::Unsigned) where {T, X <: FixedPoint{T}} = copysign(x, signed(y))
    flipsign(x::X, y::Unsigned) where {T, X <: FixedPoint{T}} = flipsign(x, signed(y))
end
@noinline function throw_not_a_signed_number_error(x)
    throw(ArgumentError("$x is not a signed number."))
end

signbit(x::X) where {X <: FixedPoint} = signbit(x.i)

~(x::X) where {X <: FixedPoint} = X(~x.i, 0)

function round(x::FixedPoint, r::RoundingMode=RoundNearest;
               digits::Union{Nothing, Integer}=nothing,
               sigdigits::Union{Nothing, Integer}=nothing,
               base::Union{Nothing, Integer}=nothing)

    d = digits === nothing ? 0 : Int(digits)

    if base !== nothing && base != 10
        throw(ArgumentError("`base` numbers other than 10 are not supported."))
    elseif sigdigits !== nothing
        throw(ArgumentError("`sigdigits` is not supported."))
    elseif d < 0
        throw(ArgumentError("negative `digits` is not supported."))
    end
    d === 0 ? _round_digits0(x, r) : _round_digits(x, r, d)
end
function _round_digits(x::X, r::RoundingMode, d::Int) where {T, f, X <: FixedPoint{T,f}}
    log10_2 = 0.3010299956639812
    d > floor(Int, ((f + 1) * log10_2)) && return x
    r = round(float(x), r, digits=d)
    typemin(X) - eps(X)/2 <= r < typemax(X) + eps(X)/2 || throw_converterror(X, r)
    clamp(r, X)
end

trunc(x::X) where {X <: FixedPoint{<:Unsigned}} = floor(x)
trunc(::Type{Ti}, x::X) where {X <: FixedPoint{<:Unsigned}, Ti <: Integer} = floor(Ti, x)

for f in (:zero, :oneunit, :one, :eps, :rawone, :rawtype, :floattype)
    @eval begin
        $f(x::FixedPoint) = $f(typeof(x))
    end
end
for f in (:(==), :<, :<=, :fld1)
    @eval begin
        $f(x::X, y::X) where {X <: FixedPoint} = $f(x.i, y.i)
    end
end
for f in (:mod1, :min, :max)
    @eval begin
        $f(x::X, y::X) where {X <: FixedPoint} = X($f(x.i, y.i), 0)
    end
end
for (m, f) in ((:(:Nearest), :round),
               (:(:ToZero), :trunc),
               (:(:Up), :ceil),
               (:(:Down), :floor))
    @eval begin
        _round_digits0(x::FixedPoint, ::RoundingMode{$m}) = $f(x)
        round(::Type{Ti}, x::FixedPoint, ::RoundingMode{$m}) where {Ti <: Integer} = $f(Ti, x)
    end
end

function length(r::StepRange{X,X}) where {X <: FixedPoint{<:ShorterThanInt}}
    start, step, stop = Int(reinterpret(r.start)), Int(reinterpret(r.step)), Int(reinterpret(r.stop))
    return div((stop - start) + step, step)
end
function length(r::StepRange{X,X}) where {X <: FixedPoint}
    start, step, stop = reinterpret(r.start), reinterpret(r.step), reinterpret(r.stop)
    return checked_div(checked_add(checked_sub(stop, start), step), step)
end
function length(r::StepRange{<:FixedPoint})
    start, step, stop = float(r.start), r.step, float(r.stop)
    return div((stop - start) + step, step)
end

hasalias(::Type) = false
hasalias(::Type{X}) where {T<:NotBiggerThanInt64, f, X<:FixedPoint{T,f}} = f isa Int

# `alias_symbol` is used to define type aliases, so we need this before we
# include "src/fixed.jl" / "src/normed.jl".
function alias_symbol(@nospecialize(X))
    Symbol(type_prefix(X), nbitsint(X), 'f', nbitsfrac(X))
end

function _alias_symbol(::Type{X}) where {X <: FixedPoint}
    if @generated
        return QuoteNode(alias_symbol(X))
    else
        return alias_symbol(X)
    end
end

function print(io::IO, x::FixedPoint{T,f}) where {T,f}
    compact = get(io, :compact, false)::Bool
    log10_2 = 0.3010299956639812
    digits = min(ceil(Int, f * log10_2), compact ? 6 : typemax(Int))
    val = round(convert(Float64, x), digits=digits)
    print(io, val)
end

@inline function showtype(io::IO, ::Type{X}) where {X <: FixedPoint}
    if hasalias(X)
        write(io, _alias_symbol(X))
    else
        print(io, X)
    end
    return nothing
end

function show(io::IO, x::FixedPoint{T,f}) where {T,f}
    compact = get(io, :compact, false)::Bool
    if compact || get(io, :typeinfo, Any) === typeof(x)
        print(io, x)
    elseif hasalias(typeof(x))
        print(io, x)
        showtype(io, typeof(x))
    else
        print(io, typeof(x), '(', x, ')')
    end
    return nothing
end

if VERSION < v"1.6.0-DEV.356" # JuliaLang/julia#36107
    function Base.showarg(io::IO, a::Array{X}, toplevel) where {X<:FixedPoint}
        toplevel || print(io, "::")
        print(io, "Array{")
        showtype(io, X)
        print(io, ",$(ndims(a))}")
        toplevel && hasalias(X) && print(io, " with eltype ", X)
    end
end

include("fixed.jl")
include("normed.jl")
include("deprecations.jl")
const UF = (N0f8, N6f10, N4f12, N2f14, N0f16)

include("arithmetic/arithmetic.jl")
using .FixedPointArithmetic
# re-export
for name in names(FixedPointArithmetic.Wrapping)
    startswith(string(name), "wrapping") || continue
    @eval export $name
end
for name in names(FixedPointArithmetic.Saturating)
    startswith(string(name), "saturating") || continue
    @eval export $name
end
for name in names(FixedPointArithmetic.Checked)
    startswith(string(name), "checked") || continue
    @eval export $name
end

# construction using the (approximate) intended value, i.e., N0f8
*(x::Real, ::Type{X}) where {X <: FixedPoint} = _convert(X, x)
Wrapping.wrapping_mul(x::Real, ::Type{X}) where {X <: FixedPoint} = x % X
Saturating.saturating_mul(x::Real, ::Type{X}) where {X <: FixedPoint} = clamp(x, X)
Checked.checked_mul(x::Real, ::Type{X}) where {X <: FixedPoint} = _convert(X, x)

# type modulus
rem(x::Real, ::Type{X}) where {X <: FixedPoint} = _rem(x, X)
Wrapping.wrapping_rem(x::Real, ::Type{X}) where {X<:FixedPoint} = _rem(x, X)
Saturating.saturating_rem(x::Real, ::Type{X}) where {X<:FixedPoint} = _rem(x, X)
Checked.checked_rem(x::Real, ::Type{X}) where {X<:FixedPoint} = _rem(x, X)

# Promotions
promote_rule(::Type{X}, ::Type{Tf}) where {X <: FixedPoint, Tf <: AbstractFloat} =
    promote_type(floattype(X), Tf)

# Note that `Tr` does not always have enough domains.
promote_rule(::Type{X}, ::Type{Tr}) where {X <: FixedPoint, Tr <: Rational} = Tr

promote_rule(::Type{X}, ::Type{Ti}) where {X <: FixedPoint, Ti <: Integer} = floattype(X)

function promote_rule(::Type{X1}, ::Type{X2}) where {T1, f1, X1 <: FixedPoint{T1,f1},
                                                     T2, f2, X2 <: FixedPoint{T2,f2}}
    X = wrapper(X1)
    X !== wrapper(X2) && return promote_type(floattype(X1), floattype(X2))

    f = max(f1, f2)  # ensure we have enough precision
    Tp = promote_type(T1, T2)
    T = (T1 <: Signed || T2 <: Signed) ? signedtype(Tp) : Tp
    # make sure we have enough integer bits
    m = max(nbitsint(X1), nbitsint(X2))
    _widen_rawtype(X{T,f}, m)
end
# TODO: avoid using @pure
@pure function _widen_rawtype(::Type{X}, m) where {T, f, X<:FixedPoint{T,f}}
    nbitsint(X) >= m && return X
    Tw = widen1(T)
    T === Tw && return X
    _widen_rawtype(wrapper(X){Tw,f}, m)
end

# Promotions for reductions
const Treduce = Float64
Base.add_sum(x::FixedPoint, y::FixedPoint) = Treduce(x) + Treduce(y)
Base.reduce_empty(::typeof(Base.add_sum), ::Type{F}) where {F<:FixedPoint}  = zero(Treduce)
Base.reduce_first(::typeof(Base.add_sum), x::FixedPoint)   = Treduce(x)
Base.mul_prod(x::FixedPoint, y::FixedPoint) = Treduce(x) * Treduce(y)
Base.reduce_empty(::typeof(Base.mul_prod), ::Type{F}) where {F<:FixedPoint} = one(Treduce)
Base.reduce_first(::typeof(Base.mul_prod), x::FixedPoint)  = Treduce(x)

"""
    sd, ad = scaledual(s::Number, a)

Return `sd` and `ad` such that `sd * ad == s * a`.
When `a` is an array of FixedPoint numbers, `sd*ad` might be faster to compute than `s*a`.
"""
scaledual(b::Number, x::Union{Number,AbstractArray{<:Number}}) = b, x
scaledual(b::Number, x::FixedPoint) = b/rawone(x), reinterpret(x)
scaledual(b::Number, x::AbstractArray{T}) where T <: FixedPoint =
    b/rawone(T), reinterpret(rawtype(T), x)

scaledual(::Type{Tdual}, x::Union{Number,AbstractArray{<:Number}}) where Tdual = oneunit(Tdual), x
scaledual(::Type{Tdual}, x::FixedPoint) where Tdual = convert(Tdual, 1/rawone(x)), reinterpret(x)
scaledual(::Type{Tdual}, x::AbstractArray{T}) where {Tdual, T <: FixedPoint} =
    convert(Tdual, 1/rawone(T)), reinterpret(rawtype(T), x)

@noinline function throw_converterror(::Type{X}, @nospecialize(x)) where X <: FixedPoint
    nbits = bitwidth(rawtype(X))
    io = IOBuffer()
    showtype(io, X)
    print(io, " is ")
    print(io, nbits == 8 ? "an " : "a ", nbits, "-bit type representing ")
    print(io, nbits <= 16 ? string(2^nbits) : "2^$nbits", " values from ")
    print(IOContext(io, :compact=>true), typemin(X), " to ")
    print(IOContext(io, :compact=>true), typemax(X), "; ")
    print(io, "cannot represent ", x)
    throw(ArgumentError(String(take!(io))))
end


function Random.rand(r::AbstractRNG, ::SamplerType{X}) where X <: FixedPoint
    X(rand(r, rawtype(X)), 0)
end

function rand!(r::AbstractRNG, A::Array{X}, ::SamplerType{X}) where {T, X <: FixedPoint{T}}
    At = unsafe_wrap(Array, reinterpret(Ptr{T}, pointer(A)), size(A))
    Random.rand!(r, At, SamplerType{T}())
    A
end

if !isdefined(Base, :get_extension)
    include("../ext/FixedPointNumbersStatisticsExt.jl")
end

if VERSION >= v"1.1" # work around https://github.com/JuliaLang/julia/issues/34121
    include("precompile.jl")
    _precompile_()
end

end # module
