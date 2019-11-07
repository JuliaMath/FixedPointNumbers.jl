# Normed{T,f} maps UInts from 0 to 2^f-1 to the range [0.0, 1.0]
# For example, Normed{UInt8,8} == N0f8 maps 0x00 to 0.0 and 0xff to 1.0

struct Normed{T<:Unsigned,f} <: FixedPoint{T,f}
    i::T

    Normed{T, f}(i::Integer,_) where {T,f} = new{T, f}(i%T)   # for setting by raw representation
end

Normed{T, f}(x::AbstractChar) where {T,f} = throw(ArgumentError("Normed cannot be constructed from a Char"))
Normed{T, f}(x::Complex) where {T,f} = Normed{T, f}(convert(real(typeof(x)), x))
Normed{T, f}(x::Base.TwicePrecision) where {T,f} = Normed{T, f}(convert(Float64, x))
Normed{T1,f}(x::Normed{T2,f}) where {T1 <: Unsigned,T2 <: Unsigned,f} = Normed{T1,f}(convert(T1, x.i), 0)

typechar(::Type{X}) where {X <: Normed} = 'N'
signbits(::Type{X}) where {X <: Normed} = 0

for T in (UInt8, UInt16, UInt32, UInt64)
    for f in 1:sizeof(T)*8
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
function Normed{T,f}(x::Normed{T2}) where {T <: Unsigned,T2 <: Unsigned,f}
    U = Normed{T,f}
    y = round((rawone(U)/rawone(x))*reinterpret(x))
    (0 <= y) & (y <= typemax(T)) || throw_converterror(U, x)
    reinterpret(U, _unsafe_trunc(T, y))
end
N0f16(x::N0f8) = reinterpret(N0f16, convert(UInt16, 0x0101*reinterpret(x)))

(::Type{U})(x::Real) where {U <: Normed} = _convert(U, x)

function _convert(::Type{U}, x) where {T, f, U <: Normed{T,f}}
    if T == UInt128 # for UInt128, we can't widen
        # the upper limit is not exact
        (0 <= x) & (x <= (typemax(T)/rawone(U))) || throw_converterror(U, x)
        y = round(rawone(U)*x)
    else
        y = round(widen1(rawone(U))*x)
        (0 <= y) & (y <= typemax(T)) || throw_converterror(U, x)
    end
    reinterpret(U, _unsafe_trunc(T, y))
end
# Prevent overflow (https://discourse.julialang.org/t/saving-greater-than-8-bit-images/6057)
function _convert(::Type{U}, x::Float16) where {T, f, U <: Normed{T,f}}
    if Float16(typemax(T)/rawone(U)) > Float32(typemax(T)/rawone(U))
        x == Float16(typemax(T)/rawone(U)) && return typemax(U)
    end
    return _convert(U, Float32(x))
end
function _convert(::Type{U}, x::Tf) where {T, f, U <: Normed{T,f}, Tf <: Union{Float32, Float64}}
    if T == UInt128 && f == 53
        0 <= x <= Tf(3.777893186295717e22) || throw_converterror(U, x)
    else
        0 <= x <= Tf((typemax(T)-rawone(U))/rawone(U)+1) || throw_converterror(U, x)
    end

    significand_bits = Tf == Float64 ? 52 : 23
    if f <= (significand_bits + 1) && sizeof(T) * 8 < significand_bits
        return reinterpret(U, unsafe_trunc(T, round(rawone(U) * x)))
    end
    # cf. the implementation of `frexp`
    Tw = f < sizeof(T) * 8 ? T : widen1(T)
    bits = sizeof(Tw) * 8 - 1
    xu = reinterpret(Tf == Float64 ? UInt64 : UInt32, x)
    k = Int(xu >> significand_bits)
    k == 0 && return zero(U) # neglect subnormal numbers
    significand = xu | (one(xu) << significand_bits)
    yh = unsafe_trunc(Tw, significand) << (bits - significand_bits)
    exponent_bias = Tf == Float64 ? 1023 : 127
    ex = exponent_bias - k + bits - f
    yi = bits >= f ? yh - (yh >> f) : yh
    if ex <= 0
        ex == 0 && return reinterpret(U, unsafe_trunc(T, yi))
        ex != -1 || signbit(signed(yi)) && return typemax(U)
        return reinterpret(U, unsafe_trunc(T, yi + yi))
    end
    ex > bits && return reinterpret(U, ex == bits + 1 ? one(T) : zero(T))
    yi += one(Tw)<<((ex - 1) & bits) # RoundNearestTiesUp
    return reinterpret(U, unsafe_trunc(T, yi >> (ex & bits)))
end

rem(x::T, ::Type{T}) where {T <: Normed} = x
rem(x::Normed, ::Type{T}) where {T <: Normed} = reinterpret(T, _unsafe_trunc(rawtype(T), round((rawone(T)/rawone(x))*reinterpret(x))))
rem(x::Real, ::Type{T}) where {T <: Normed} = reinterpret(T, _unsafe_trunc(rawtype(T), round(rawone(T)*x)))
rem(x::Float16, ::Type{T}) where {T <: Normed} = rem(Float32(x), T)  # avoid overflow

float(x::Normed) = convert(floattype(x), x)

Base.BigFloat(x::Normed) = reinterpret(x)*(1/BigFloat(rawone(x)))
function (::Type{T})(x::Normed) where {T <: AbstractFloat}
    y = reinterpret(x)*(one(rawtype(x))/convert(T, rawone(x)))
    convert(T, y)  # needed for types like Float16 which promote arithmetic to Float32
end
Base.Bool(x::Normed) = x == zero(x) ? false : true
Base.Integer(x::Normed) = convert(Integer, x*1.0)
(::Type{T})(x::Normed) where {T <: Integer} = convert(T, x*(1/oneunit(T)))
Base.Rational{Ti}(x::Normed) where {Ti <: Integer} = convert(Ti, reinterpret(x))//convert(Ti, rawone(x))
Base.Rational(x::Normed) = reinterpret(x)//rawone(x)

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
# The main subtlety here is that iterating over N0f8(0):N0f8(1) will wrap around
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
        Tw = widen1(T)
        T == Tw && break
        T = Tw
        i = 8*sizeof(T)-f
    end
    :(Normed{$T,$f})
end

_unsafe_trunc(::Type{T}, x::Integer) where {T} = x % T
_unsafe_trunc(::Type{T}, x) where {T}          = unsafe_trunc(T, x)
if !signbit(signed(unsafe_trunc(UInt, -12.345)))
    # a workaround for 32-bit ARMv7 (issue #134)
    function _unsafe_trunc(::Type{T}, x::AbstractFloat) where {T}
        unsafe_trunc(T, unsafe_trunc(typeof(signed(zero(T))), x))
    end
end
