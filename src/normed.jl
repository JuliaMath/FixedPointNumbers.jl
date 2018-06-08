# Normed{T,f} maps UInts from 0 to 2^f-1 to the range [0.0, 1.0]
# For example, Normed{UInt8,8} == N0f8 maps 0x00 to 0.0 and 0xff to 1.0

struct Normed{T<:Unsigned,f} <: FixedPoint{T,f}
    i::T

    Normed{T, f}(i::Integer,_) where {T,f} = new{T, f}(i%T)   # for setting by raw representation
    Normed{T, f}(x) where {T,f} = convert(Normed{T,f}, x)
    Normed{T, f}(x::Normed{T,f}) where {T,f} = x
    Normed{T, f}(x::AbstractChar) where {T,f} = throw(ArgumentError("Normed cannot be constructed from a Char"))
    Normed{T, f}(x::Complex) where {T,f} = Normed{T, f}(convert(real(typeof(x)), x))
    Normed{T, f}(x::Base.TwicePrecision) where {T,f} = Normed{T, f}(convert(Float64, x))
end

typechar(::Type{X}) where {X <: Normed} = 'N'
signbits(::Type{X}) where {X <: Normed} = 0

for T in (UInt8, UInt16, UInt32, UInt64)
    for f in 0:sizeof(T)*8
        sym = Symbol(String(take!(showtype(_iotypealias, Normed{T,f}))))
        @eval begin
            const $sym = Normed{$T,$f}
            export $sym
        end
    end
end

reinterpret(::Type{Normed{T,f}}, x::T) where {T <: Unsigned,f} = Normed{T,f}(x, 0)

zero(::Type{Normed{T,f}}) where {T,f} = Normed{T,f}(zero(T),0)
function oneunit(::Type{T}) where {T <: Normed}
    T(typemax(rawtype(T)) >> (8*sizeof(T)-nbitsfrac(T)), 0)
end
one(::Type{T}) where {T <: Normed} = oneunit(T)
zero(x::Normed) = zero(typeof(x))
oneunit(x::Normed) =  one(typeof(x))
one(x::Normed) = oneunit(x)
rawone(v) = reinterpret(one(v))

# Conversions
convert(::Type{U}, x::U) where {U <: Normed} = x
convert(::Type{Normed{T1,f}}, x::Normed{T2,f}) where {T1 <: Unsigned,T2 <: Unsigned,f} = Normed{T1,f}(convert(T1, x.i), 0)
function convert(::Type{Normed{T,f}}, x::Normed{T2}) where {T <: Unsigned,T2 <: Unsigned,f}
    U = Normed{T,f}
    y = round((rawone(U)/rawone(x))*reinterpret(x))
    (0 <= y) & (y <= typemax(T)) || throw_converterror(U, x)
    reinterpret(U, _unsafe_trunc(T, y))
end
convert(::Type{U}, x::Real) where {U <: Normed} = _convert(U, rawtype(U), x)

convert(::Type{N0f16}, x::N0f8) = reinterpret(N0f16, convert(UInt16, 0x0101*reinterpret(x)))
function _convert(::Type{U}, ::Type{T}, x) where {U <: Normed,T}
    y = round(widen1(rawone(U))*x)
    (0 <= y) & (y <= typemax(T)) || throw_converterror(U, x)
    U(_unsafe_trunc(T, y), 0)
end
# Prevent overflow (https://discourse.julialang.org/t/saving-greater-than-8-bit-images/6057)
_convert(::Type{U}, ::Type{T}, x::Float16) where {U <: Normed,T} =
    _convert(U, T, Float32(x))
_convert(::Type{U}, ::Type{UInt128}, x::Float16) where {U <: Normed} =
    _convert(U, UInt128, Float32(x))
function _convert(::Type{U}, ::Type{UInt128}, x) where {U <: Normed}
    y = round(rawone(U)*x)   # for UInt128, we can't widen
    (0 <= y) & (y <= typemax(UInt128)) & (x <= convert(Float64, typemax(U))) || throw_converterror(U, x)
    U(_unsafe_trunc(UInt128, y), 0)
end

rem(x::T, ::Type{T}) where {T <: Normed} = x
rem(x::Normed, ::Type{T}) where {T <: Normed} = reinterpret(T, _unsafe_trunc(rawtype(T), round((rawone(T)/rawone(x))*reinterpret(x))))
rem(x::Real, ::Type{T}) where {T <: Normed} = reinterpret(T, _unsafe_trunc(rawtype(T), round(rawone(T)*x)))
rem(x::Float16, ::Type{T}) where {T <: Normed} = rem(Float32(x), T)  # avoid overflow

# convert(::Type{AbstractFloat}, x::Normed) = convert(floattype(x), x)
float(x::Normed) = convert(floattype(x), x)

convert(::Type{BigFloat}, x::Normed) = reinterpret(x)*(1/BigFloat(rawone(x)))
function convert(::Type{T}, x::Normed) where {T <: AbstractFloat}
    y = reinterpret(x)*(one(rawtype(x))/convert(T, rawone(x)))
    convert(T, y)  # needed for types like Float16 which promote arithmetic to Float32
end
convert(::Type{Bool}, x::Normed) = x == zero(x) ? false : true
convert(::Type{Integer}, x::Normed) = convert(Integer, x*1.0)
convert(::Type{T}, x::Normed) where {T <: Integer} = convert(T, x*(1/oneunit(T)))
convert(::Type{Rational{Ti}}, x::Normed) where {Ti <: Integer} = convert(Ti, reinterpret(x))//convert(Ti, rawone(x))
convert(::Type{Rational}, x::Normed) = reinterpret(x)//rawone(x)

# Traits
abs(x::Normed) = x

(-)(x::T) where {T <: Normed} = T(-reinterpret(x), 0)
(~)(x::T) where {T <: Normed} = T(~reinterpret(x), 0)

+(x::Normed{T,f}, y::Normed{T,f}) where {T,f} = Normed{T,f}(convert(T, x.i+y.i),0)
-(x::Normed{T,f}, y::Normed{T,f}) where {T,f} = Normed{T,f}(convert(T, x.i-y.i),0)
*(x::T, y::T) where {T <: Normed} = convert(T,convert(floattype(T), x)*convert(floattype(T), y))
/(x::T, y::T) where {T <: Normed} = convert(T,convert(floattype(T), x)/convert(floattype(T), y))

# Comparisons
 <(x::T, y::T) where {T <: Normed} = reinterpret(x) < reinterpret(y)
<=(x::T, y::T) where {T <: Normed} = reinterpret(x) <= reinterpret(y)

# Functions
trunc(x::T) where {T <: Normed} = T(div(reinterpret(x), rawone(T))*rawone(T),0)
floor(x::T) where {T <: Normed} = trunc(x)
function round(x::Normed{T,f}) where {T,f}
    mask = convert(T, 1<<(f-1))
    y = trunc(x)
    return convert(T, reinterpret(x)-reinterpret(y)) & mask>0 ?
            Normed{T,f}(y+oneunit(Normed{T,f})) : y
end
function ceil(x::Normed{T,f}) where {T,f}
    k = 8*sizeof(T)-f
    mask = (typemax(T)<<k)>>k
    y = trunc(x)
    return convert(T, reinterpret(x)-reinterpret(y)) & (mask)>0 ?
            Normed{T,f}(y+oneunit(Normed{T,f})) : y
end

trunc(::Type{T}, x::Normed) where {T <: Integer} = convert(T, div(reinterpret(x), rawone(x)))
round(::Type{T}, x::Normed) where {T <: Integer} = round(T, reinterpret(x)/rawone(x))
floor(::Type{T}, x::Normed) where {T <: Integer} = trunc(T, x)
 ceil(::Type{T}, x::Normed) where {T <: Integer} =  ceil(T, reinterpret(x)/rawone(x))

isfinite(x::Normed) = true
isnan(x::Normed) = false
isinf(x::Normed) = false

bswap(x::Normed{UInt8,f}) where {f} = x
bswap(x::Normed)  = typeof(x)(bswap(reinterpret(x)),0)

function minmax(x::T, y::T) where {T <: Normed}
    a, b = minmax(reinterpret(x), reinterpret(y))
    T(a,0), T(b,0)
end

# Iteration
# The main subtlety here is that iterating over 0x00uf8:0xffuf8 will wrap around
# unless we iterate using a wider type
@inline start(r::StepRange{T}) where {T <: Normed} = widen1(reinterpret(r.start))
@inline next(r::StepRange{T}, i::Integer) where {T <: Normed} = (T(i,0), i+reinterpret(r.step))
@inline function done(r::StepRange{T}, i::Integer) where {T <: Normed}
    i1, i2 = reinterpret(r.start), reinterpret(r.stop)
    isempty(r) | (i < min(i1, i2)) | (i > max(i1, i2))
end

function decompose(x::Normed)
    g = gcd(reinterpret(x), rawone(x))
    div(reinterpret(x),g), 0, div(rawone(x),g)
end

# Promotions
promote_rule(::Type{T}, ::Type{Tf}) where {T <: Normed,Tf <: AbstractFloat} = promote_type(floattype(T), Tf)
promote_rule(::Type{T}, ::Type{R}) where {T <: Normed,R <: Rational} = R
function promote_rule(::Type{T}, ::Type{Ti}) where {T <: Normed,Ti <: Union{Signed, Unsigned}}
    floattype(T)
end
@generated function promote_rule(::Type{Normed{T1,f1}}, ::Type{Normed{T2,f2}}) where {T1,T2,f1,f2}
    f = max(f1, f2)  # ensure we have enough precision
    T = promote_type(T1, T2)
    # make sure we have enough integer bits
    i1, i2 = 8*sizeof(T1)-f1, 8*sizeof(T2)-f2  # number of integer bits for each
    i = 8*sizeof(T)-f
    while i < max(i1, i2)
        T = widen1(T)
        i = 8*sizeof(T)-f
    end
    :(Normed{$T,$f})
end

_unsafe_trunc(::Type{T}, x::Integer) where {T} = x % T
_unsafe_trunc(::Type{T}, x) where {T}          = unsafe_trunc(T, x)
