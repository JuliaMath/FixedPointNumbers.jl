"""
Module for utility functions and macros, which are independent of `FixedPoint`.
This also includes the declaration of public or internal APIs.
"""
module Utilities

export ShortInts
export LongInts
export ShorterThanInt
export NotBiggerThanInt
export NotBiggerThanInt64
export SShorterThanInt
export UShorterThanInt

export bitwidth, widen1, signedtype, wrapper
export @f32, @exp2, significand_bits, exponent_bias
export _unsafe_trunc
export div_2f, div_2fm1

export floattype, rawone, inv_rawone, nbitsfrac, rawtype, signbits, nbitsint
export scaledual


const ShortInts = Union{Int8, UInt8, Int16, UInt16}
const LongInts = Union{Int64, UInt64, Int128, UInt128, BigInt}

const ShorterThanInt = Int === Int32 ? ShortInts : Union{ShortInts, Int32, UInt32}
const NotBiggerThanInt = Union{ShorterThanInt, Int, UInt}
const NotBiggerThanInt64 = Union{ShortInts, Int32, UInt32, Int64, UInt64}
const SShorterThanInt = typeintersect(ShorterThanInt, Signed)
const UShorterThanInt = typeintersect(ShorterThanInt, Unsigned)

bitwidth(T::Type) = 8sizeof(T)

widen1(T::Type)         = T # fallback
widen1(::Type{Int8})    = Int16
widen1(::Type{UInt8})   = UInt16
widen1(::Type{Int16})   = Int32
widen1(::Type{UInt16})  = UInt32
widen1(::Type{Int32})   = Int64
widen1(::Type{UInt32})  = UInt64
widen1(::Type{Int64})   = Int128
widen1(::Type{UInt64})  = UInt128
widen1(x::Integer) = x % widen1(typeof(x))

signedtype(::Type{T}) where {T <: Integer} = typeof(signed(zero(T)))

wrapper(@nospecialize(T)) = Base.typename(T).wrapper

macro f32(x::Float64) # just for hexadecimal floating-point literals
    :(Float32($x))
end
macro exp2(n)
     :(_exp2(Val($(esc(n)))))
end
_exp2(::Val{N}) where {N} = exp2(N)

# these are defined in julia/float.jl or julia/math.jl, but not exported
significand_bits(::Type{Float32}) = 23
significand_bits(::Type{Float64}) = 52
exponent_bias(::Type{Float32}) = 127
exponent_bias(::Type{Float64}) = 1023

_unsafe_trunc(::Type{T}, x::Integer) where {T} = x % T
_unsafe_trunc(::Type{T}, x) where {T} = unsafe_trunc(T, x)
# issue #202, #211
_unsafe_trunc(::Type{T}, x::BigFloat) where {T <: Integer} = trunc(BigInt, x) % T

# issue #288
function _unsafe_trunc(::Type{T}, x::AbstractFloat) where {T <: Integer}
    if T <: ShortInts
        return unsafe_trunc(Int32, x) % T
    elseif T <: Unsigned
        return copysign(unsafe_trunc(T, abs(x)), x)
    else
        return unsafe_trunc(T, x)
    end
end

# Division by `2^f` with RoundNearest.
function div_2f(x::T, ::Val{f}) where {T,f}
    xf = x & (T(-1) >>> (bitwidth(T) - f - 1))
    half = oneunit(T) << (f - 1)
    c = half - (xf === half)
    (x + c) >> f
end
div_2f(x::T, ::Val{0}) where {T} = x

# Division by `2^f-1` with RoundNearest. The result would be in the lower half bits.
div_2fm1(x::T, ::Val{f}) where {T,f} = (x + (T(1) << (f - 1) - 0x1)) ÷ (T(1) << f - 0x1)
div_2fm1(x::T, ::Val{1}) where {T} = x
div_2fm1(x::UInt16, ::Val{8}) = (((x + 0x80) >> 0x8) + x + 0x80) >> 0x8
div_2fm1(x::UInt32, ::Val{16}) = (((x + 0x8000) >> 0x10) + x + 0x8000) >> 0x10
div_2fm1(x::UInt64, ::Val{32}) = (((x + 0x80000000) >> 0x20) + x + 0x80000000) >> 0x20
div_2fm1(x::UInt128, ::Val{64}) = (((x + 0x8000000000000000) >> 0x40) + x + 0x8000000000000000) >> 0x40



floattype(::Type{T}) where {T<:AbstractFloat} = T # fallback (we want a MethodError if no method producing AbstractFloat is defined)
floattype(::Type{T}) where {T<:Union{ShortInts,Bool}} = Float32
floattype(::Type{T}) where {T<:Integer} = Float64
floattype(::Type{T}) where {T<:LongInts} = BigFloat
floattype(::Type{T}) where {I<:Integer,T<:Rational{I}} = typeof(zero(I) / oneunit(I))
floattype(::Type{<:AbstractIrrational}) = Float64
# Non-Real types
floattype(::Type{Complex{T}}) where {T} = Complex{floattype(T)}
floattype(::Type{Base.TwicePrecision{Float64}}) = Float64    # wider would be nice, but hardware support is paramount
floattype(::Type{Base.TwicePrecision{T}}) where {T<:Union{Float16,Float32}} = widen(T)

function rawone end

# for Julia v1.0, which does not fold `div_float` before inlining
inv_rawone(x) = (@generated) ? (y = 1.0 / rawone(x); :($y)) : 1.0 / rawone(x)

function nbitsfrac end
function rawtype end
function signbits end
function nbitsint end

function scaledual end

end # module Utilities
