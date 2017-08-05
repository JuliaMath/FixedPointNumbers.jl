# Normed{T,f} maps UInts from 0 to 2^f-1 to the range [0.0, 1.0]
# For example, Normed{UInt8,8} == N0f8 maps 0x00 to 0.0 and 0xff to 1.0

struct Normed{T<:Unsigned,f} <: FixedPoint{T,f}
    i::T

    (::Type{Normed{T, f}}){T, f}(i::Integer,_) = new{T, f}(i%T)   # for setting by raw representation
    (::Type{Normed{T, f}}){T, f}(x) = convert(Normed{T,f}, x)
end

typechar{X<:Normed}(::Type{X}) = 'N'
signbits{X<:Normed}(::Type{X}) = 0

for T in (UInt8, UInt16, UInt32, UInt64)
    for f in 0:sizeof(T)*8
        sym = Symbol(String(take!(showtype(_iotypealias, Normed{T,f}))))
        @eval begin
            const $sym = Normed{$T,$f}
            export $sym
        end
    end
end

reinterpret{T<:Unsigned, f}(::Type{Normed{T,f}}, x::T) = Normed{T,f}(x, 0)

zero{T,f}(::Type{Normed{T,f}}) = Normed{T,f}(zero(T),0)
function one{T<:Normed}(::Type{T})
    T(typemax(rawtype(T)) >> (8*sizeof(T)-nbitsfrac(T)), 0)
end
zero(x::Normed) = zero(typeof(x))
 one(x::Normed) =  one(typeof(x))
rawone(v) = reinterpret(one(v))

# Conversions
convert{U<:Normed}(::Type{U}, x::U) = x
convert{T1<:Unsigned,T2<:Unsigned,f}(::Type{Normed{T1,f}}, x::Normed{T2,f}) = Normed{T1,f}(convert(T1, x.i), 0)
function convert{T<:Unsigned,T2<:Unsigned,f}(::Type{Normed{T,f}}, x::Normed{T2})
    U = Normed{T,f}
    y = round((rawone(U)/rawone(x))*reinterpret(x))
    (0 <= y) & (y <= typemax(T)) || throw_converterror(U, x)
    reinterpret(U, _unsafe_trunc(T, y))
end
convert{U<:Normed}(::Type{U}, x::Real) = _convert(U, rawtype(U), x)

convert(::Type{N0f16}, x::N0f8) = reinterpret(N0f16, convert(UInt16, 0x0101*reinterpret(x)))
function _convert{U<:Normed,T}(::Type{U}, ::Type{T}, x)
    y = round(widen1(rawone(U))*x)
    (0 <= y) & (y <= typemax(T)) || throw_converterror(U, x)
    U(_unsafe_trunc(T, y), 0)
end
function _convert{U<:Normed}(::Type{U}, ::Type{UInt128}, x)
    y = round(rawone(U)*x)   # for UInt128, we can't widen
    (0 <= y) & (y <= typemax(UInt128)) & (x <= Float64(typemax(U))) || throw_converterror(U, x)
    U(_unsafe_trunc(UInt128, y), 0)
end

rem{T<:Normed}(x::T, ::Type{T}) = x
rem{T<:Normed}(x::Normed, ::Type{T}) = reinterpret(T, _unsafe_trunc(rawtype(T), round((rawone(T)/rawone(x))*reinterpret(x))))
rem{T<:Normed}(x::Real, ::Type{T}) = reinterpret(T, _unsafe_trunc(rawtype(T), round(rawone(T)*x)))

# convert(::Type{AbstractFloat}, x::Normed) = convert(floattype(x), x)
float(x::Normed) = convert(floattype(x), x)

convert(::Type{BigFloat}, x::Normed) = reinterpret(x)*(1/BigFloat(rawone(x)))
function convert{T<:AbstractFloat}(::Type{T}, x::Normed)
    y = reinterpret(x)*(one(rawtype(x))/convert(T, rawone(x)))
    convert(T, y)  # needed for types like Float16 which promote arithmetic to Float32
end
convert(::Type{Bool}, x::Normed) = x == zero(x) ? false : true
convert(::Type{Integer}, x::Normed) = convert(Integer, x*1.0)
convert{T<:Integer}(::Type{T}, x::Normed) = convert(T, x*(1/one(T)))
convert{Ti<:Integer}(::Type{Rational{Ti}}, x::Normed) = convert(Ti, reinterpret(x))//convert(Ti, rawone(x))
convert(::Type{Rational}, x::Normed) = reinterpret(x)//rawone(x)

# Traits
abs(x::Normed) = x

(-){T<:Normed}(x::T) = T(-reinterpret(x), 0)
(~){T<:Normed}(x::T) = T(~reinterpret(x), 0)

+{T,f}(x::Normed{T,f}, y::Normed{T,f}) = Normed{T,f}(convert(T, x.i+y.i),0)
-{T,f}(x::Normed{T,f}, y::Normed{T,f}) = Normed{T,f}(convert(T, x.i-y.i),0)
*{T<:Normed}(x::T, y::T) = convert(T,convert(floattype(T), x)*convert(floattype(T), y))
/{T<:Normed}(x::T, y::T) = convert(T,convert(floattype(T), x)/convert(floattype(T), y))

# Comparisons
 <{T<:Normed}(x::T, y::T) = reinterpret(x) < reinterpret(y)
<={T<:Normed}(x::T, y::T) = reinterpret(x) <= reinterpret(y)

# Functions
trunc{T<:Normed}(x::T) = T(div(reinterpret(x), rawone(T))*rawone(T),0)
floor{T<:Normed}(x::T) = trunc(x)
function round{T,f}(x::Normed{T,f})
    mask = convert(T, 1<<(f-1))
    y = trunc(x)
    return convert(T, reinterpret(x)-reinterpret(y)) & mask>0 ?
            Normed{T,f}(y+one(Normed{T,f})) : y
end
function ceil{T,f}(x::Normed{T,f})
    k = 8*sizeof(T)-f
    mask = (typemax(T)<<k)>>k
    y = trunc(x)
    return convert(T, reinterpret(x)-reinterpret(y)) & (mask)>0 ?
            Normed{T,f}(y+one(Normed{T,f})) : y
end

trunc{T<:Integer}(::Type{T}, x::Normed) = convert(T, div(reinterpret(x), rawone(x)))
round{T<:Integer}(::Type{T}, x::Normed) = round(T, reinterpret(x)/rawone(x))
floor{T<:Integer}(::Type{T}, x::Normed) = trunc(T, x)
 ceil{T<:Integer}(::Type{T}, x::Normed) =  ceil(T, reinterpret(x)/rawone(x))

isfinite(x::Normed) = true
isnan(x::Normed) = false
isinf(x::Normed) = false

bswap{f}(x::Normed{UInt8,f}) = x
bswap(x::Normed)  = typeof(x)(bswap(reinterpret(x)),0)

function minmax{T<:Normed}(x::T, y::T)
    a, b = minmax(reinterpret(x), reinterpret(y))
    T(a,0), T(b,0)
end

# Iteration
# The main subtlety here is that iterating over 0x00uf8:0xffuf8 will wrap around
# unless we iterate using a wider type
@inline start{T<:Normed}(r::StepRange{T}) = widen1(reinterpret(r.start))
@inline next{T<:Normed}(r::StepRange{T}, i::Integer) = (T(i,0), i+reinterpret(r.step))
@inline function done{T<:Normed}(r::StepRange{T}, i::Integer)
    i1, i2 = reinterpret(r.start), reinterpret(r.stop)
    isempty(r) | (i < min(i1, i2)) | (i > max(i1, i2))
end

function decompose(x::Normed)
    g = gcd(reinterpret(x), rawone(x))
    div(reinterpret(x),g), 0, div(rawone(x),g)
end

# Promotions
promote_rule{T<:Normed,Tf<:AbstractFloat}(::Type{T}, ::Type{Tf}) = promote_type(floattype(T), Tf)
promote_rule{T<:Normed, R<:Rational}(::Type{T}, ::Type{R}) = R
function promote_rule{T<:Normed, Ti<:Union{Signed,Unsigned}}(::Type{T}, ::Type{Ti})
    floattype(T)
end
@generated function promote_rule{T1,T2,f1,f2}(::Type{Normed{T1,f1}}, ::Type{Normed{T2,f2}})
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

_unsafe_trunc{T}(::Type{T}, x::Integer) = x % T
_unsafe_trunc{T}(::Type{T}, x)          = unsafe_trunc(T, x)
