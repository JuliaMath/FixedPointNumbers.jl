module FixedPointNumbers

import Base: ==, <, <=, -, +, *, /, ~, isapprox,
             convert, promote_rule, show, bitstring, abs, decompose,
             isnan, isinf, isfinite, isinteger,
             zero, oneunit, one, typemin, typemax, floatmin, floatmax, eps, reinterpret,
             big, rationalize, float, trunc, round, floor, ceil, bswap, clamp,
             div, fld, cld, rem, mod, mod1, fld1, min, max, minmax,
             signed, unsigned, copysign, flipsign, signbit,
             length

import Statistics   # for _mean_promote
import Random: Random, AbstractRNG, SamplerType, rand!

import Base.Checked: checked_mul, checked_div
using Base.Checked: checked_add, checked_sub

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

# construction using the (approximate) intended value, i.e., N0f8
*(x::Real, ::Type{X}) where {X <: FixedPoint} = _convert(X, x)

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

For FixedPoint numbers, the default criterion is that `x` and `y` differ by no more than `eps`, the separation between adjacent fixed-point numbers.
"""
function isapprox(x::T, y::T; rtol=0, atol=max(eps(x), eps(y))) where {T <: FixedPoint}
    maxdiff = T(atol+rtol*max(abs(x), abs(y)))
    rx, ry, rd = reinterpret(x), reinterpret(y), reinterpret(maxdiff)
    abs(signed(widen1(rx))-signed(widen1(ry))) <= rd
end
function isapprox(x::FixedPoint, y::FixedPoint; rtol=0, atol=max(eps(x), eps(y)))
    isapprox(promote(x, y)...; rtol=rtol, atol=atol)
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

# for Julia v1.0, which does not fold `div_float` before inlining
inv_rawone(x) = (@generated) ? (y = 1.0 / rawone(x); :($y)) : 1.0 / rawone(x)

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
floattype(::Type{T}) where {T <: AbstractFloat} = T # fallback (we want a MethodError if no method producing AbstractFloat is defined)
floattype(::Type{T}) where {T <: Union{ShortInts, Bool}} = Float32
floattype(::Type{T}) where {T <: Integer} = Float64
floattype(::Type{T}) where {T <: LongInts} = BigFloat
floattype(::Type{T}) where {I <: Integer, T <: Rational{I}} = typeof(zero(I)/oneunit(I))
floattype(::Type{<:AbstractIrrational}) = Float64
floattype(::Type{X}) where {T <: ShortInts, X <: FixedPoint{T}} = Float32
floattype(::Type{X}) where {T <: Integer, X <: FixedPoint{T}} = Float64
floattype(::Type{X}) where {T <: LongInts, X <: FixedPoint{T}} = BigFloat

# Non-Real types
floattype(::Type{Complex{T}}) where T = Complex{floattype(T)}
floattype(::Type{Base.TwicePrecision{Float64}}) = Float64    # wider would be nice, but hardware support is paramount
floattype(::Type{Base.TwicePrecision{T}}) where T<:Union{Float16,Float32} = widen(T)

float(x::FixedPoint) = convert(floattype(x), x)

wrapping_mul(x::X, y::X) where {X <: FixedPoint} = (float(x) * float(y)) % X
*(x::X, y::X) where {X <: FixedPoint} = wrapping_mul(x, y)

function checked_div(x::X, y::X, r::RoundingMode = RoundToZero) where {T, X <: FixedPoint{T}}
    y === zero(X) && throw(DivideError())
    z = round(floattype(X)(x.i) / floattype(X)(y.i), r)
    if T <: Signed
        z <= typemax(T) || throw_overflowerror_div(r, x, y)
    end
    _unsafe_trunc(T, z)
end
div(x::X, y::X, r::RoundingMode = RoundToZero) where {X <: FixedPoint} = checked_div(x, y, r)
fld(x::X, y::X) where {X <: FixedPoint} = checked_div(x, y, RoundDown)
cld(x::X, y::X) where {X <: FixedPoint} = checked_div(x, y, RoundUp)

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

# Workaround for poor promotion due to lack of PR #207
function clamp(x::AbstractFloat, ::Type{X}) where {X <: FixedPoint}
    Tf = promote_type(typeof(x), floattype(X))
    clamp(Tf(x), typemin(X), typemax(X)) % X
end

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
for f in (:-, :~, :abs)
    @eval begin
        $f(x::X) where {X <: FixedPoint} = X($f(x.i), 0)
    end
end
for f in (:+, :-, :rem, :mod, :mod1, :min, :max)
    @eval begin
        $f(x::X, y::X) where {X <: FixedPoint} = X($f(x.i, y.i), 0)
    end
end
for (m, f) in ((:(:Nearest), :round),
               (:(:ToZero), :trunc),
               (:(:Up), :ceil),
               (:(:Down), :floor))
    @eval begin
        round(x::FixedPoint, ::RoundingMode{$m}) = $f(x)
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

# Printing. These are used to generate type-symbols, so we need them
# before we include any files.
function showtype(io::IO, ::Type{X}) where {X <: FixedPoint}
    if hasalias(X)
        f = nbitsfrac(X)
        m = nbitsint(X)
        write(io, typechar(X))
        m > 9 && write(io, Char(m ÷ 10 + 0x30))
        write(io, Char(m % 10 + 0x30), 'f')
        f > 9 && write(io, Char(f ÷ 10 + 0x30))
        write(io, Char(f % 10 + 0x30))
    else
        print(io, X)
    end
    io
end

function show(io::IO, x::FixedPoint{T,f}) where {T,f}
    compact = get(io, :compact, false)::Bool
    log10_2 = 0.3010299956639812
    val = round(convert(Float64, x), digits=ceil(Int, f * log10_2))
    if compact
        show(io, val)
    elseif hasalias(typeof(x))
        show(io, val)
        showtype(io, typeof(x))
    else
        print(io, typeof(x), '(', val, ')')
    end
end

function Base.showarg(io::IO, a::Array{T}, toplevel) where {T<:FixedPoint}
    toplevel || print(io, "::")
    print(io, "Array{")
    showtype(io, T)
    print(io, ",$(ndims(a))}")
    toplevel && print(io, " with eltype ", T)
end

include("fixed.jl")
include("normed.jl")
include("deprecations.jl")
const UF = (N0f8, N6f10, N4f12, N2f14, N0f16)

# Promotions for reductions
const Treduce = Float64
Base.add_sum(x::FixedPoint, y::FixedPoint) = Treduce(x) + Treduce(y)
Base.reduce_empty(::typeof(Base.add_sum), ::Type{F}) where {F<:FixedPoint}  = zero(Treduce)
Base.reduce_first(::typeof(Base.add_sum), x::FixedPoint)   = Treduce(x)
Base.mul_prod(x::FixedPoint, y::FixedPoint) = Treduce(x) * Treduce(y)
Base.reduce_empty(::typeof(Base.mul_prod), ::Type{F}) where {F<:FixedPoint} = one(Treduce)
Base.reduce_first(::typeof(Base.mul_prod), x::FixedPoint)  = Treduce(x)

if isdefined(Statistics, :_mean_promote)
    Statistics._mean_promote(x::Real, y::FixedPoint) = Treduce(y)
end

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

@noinline function throw_overflowerror(op::Symbol, @nospecialize(x), @nospecialize(y))
    io = IOBuffer()
    print(io, x, ' ', op, ' ', y, " overflowed for type ")
    showtype(io, typeof(x))
    throw(OverflowError(String(take!(io))))
end
@noinline function throw_overflowerror_div(r::RoundingMode, @nospecialize(x), @nospecialize(y))
    io = IOBuffer()
    op = r === RoundUp ? "cld(" : r === RoundDown ? "fld(" : "div("
    print(io, op, x, ", ", y, ") overflowed for type ", rawtype(x))
    throw(OverflowError(String(take!(io))))
end

function Random.rand(r::AbstractRNG, ::SamplerType{X}) where X <: FixedPoint
    X(rand(r, rawtype(X)), 0)
end

function rand!(r::AbstractRNG, A::Array{X}, ::SamplerType{X}) where {T, X <: FixedPoint{T}}
    At = unsafe_wrap(Array, reinterpret(Ptr{T}, pointer(A)), size(A))
    Random.rand!(r, At, SamplerType{T}())
    A
end

if VERSION >= v"1.1" # work around https://github.com/JuliaLang/julia/issues/34121
    include("precompile.jl")
    _precompile_()
end

end # module
