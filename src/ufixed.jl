# UFixed{T,f} maps UInts from 0 to 2^f-1 to the range [0.0, 1.0]
# For example, UFixed{UInt8,8} == N0f8 maps 0x00 to 0.0 and 0xff to 1.0

immutable UFixed{T<:Unsigned,f} <: FixedPoint{T,f}
    i::T

    UFixed(i::Integer,_) = new(i%T)   # for setting by raw representation
    UFixed(x) = convert(UFixed{T,f}, x)
end

  rawtype{T,f}(::Type{UFixed{T,f}}) = T
  rawtype(x::Number) = rawtype(typeof(x))
nbitsfrac{T,f}(::Type{UFixed{T,f}}) = f
floattype{T<:UFixed}(::Type{T}) = floattype(supertype(T))
typechar{X<:UFixed}(::Type{X}) = 'N'
signbits{X<:UFixed}(::Type{X}) = 0

for T in (UInt8, UInt16, UInt32, UInt64)
    for f in 0:sizeof(T)*8
        sym = Symbol(takebuf_string(showtype(_iotypealias, UFixed{T,f})))
        @eval begin
            typealias $sym UFixed{$T,$f}
            export $sym
        end
    end
end

reinterpret{T<:Unsigned, f}(::Type{UFixed{T,f}}, x::T) = UFixed{T,f}(x, 0)

## The next lines mimic the floating-point literal syntax "3.2f0"
# construction using a UInt, i.e., 0xccuf8
immutable UFixedConstructor{T,f} end
*{T,f}(n::Integer, ::UFixedConstructor{T,f}) = UFixed{T,f}(n,0)
const uf8  = UFixedConstructor{UInt8,8}()
const uf10 = UFixedConstructor{UInt16,10}()
const uf12 = UFixedConstructor{UInt16,12}()
const uf14 = UFixedConstructor{UInt16,14}()
const uf16 = UFixedConstructor{UInt16,16}()

zero{T,f}(::Type{UFixed{T,f}}) = UFixed{T,f}(zero(T),0)
function one{T<:UFixed}(::Type{T})
    T(typemax(rawtype(T)) >> (8*sizeof(T)-nbitsfrac(T)), 0)
end
zero(x::UFixed) = zero(typeof(x))
 one(x::UFixed) =  one(typeof(x))
rawone(v) = reinterpret(one(v))

# Conversions
convert{T<:UFixed}(::Type{T}, x::T) = x
convert{T1,T2,f}(::Type{UFixed{T1,f}}, x::UFixed{T2,f}) = UFixed{T1,f}(convert(T1, x.i), 0)
function convert{T,T2,f}(::Type{UFixed{T,f}}, x::UFixed{T2})
    U = UFixed{T,f}
    y = round((rawone(U)/rawone(x))*reinterpret(x))
    (0 <= y) & (y <= typemax(T)) || throw_converterror(U, x)
    reinterpret(U, _unsafe_trunc(T, y))
end
convert(::Type{N0f16}, x::N0f8) = reinterpret(N0f16, convert(UInt16, 0x0101*reinterpret(x)))
convert{U<:UFixed}(::Type{U}, x::Real) = _convert(U, rawtype(U), x)
function _convert{U<:UFixed,T}(::Type{U}, ::Type{T}, x)
    y = round(widen1(rawone(U))*x)
    (0 <= y) & (y <= typemax(T)) || throw_converterror(U, x)
    U(_unsafe_trunc(T, y), 0)
end
function _convert{U<:UFixed}(::Type{U}, ::Type{UInt128}, x)
    y = round(rawone(U)*x)   # for UInt128, we can't widen
    (0 <= y) & (y <= typemax(UInt128)) & (x <= Float64(typemax(U))) || throw_converterror(U, x)
    U(_unsafe_trunc(UInt128, y), 0)
end

rem{T<:UFixed}(x::T, ::Type{T}) = x
rem{T<:UFixed}(x::UFixed, ::Type{T}) = reinterpret(T, _unsafe_trunc(rawtype(T), round((rawone(T)/rawone(x))*reinterpret(x))))
rem{T<:UFixed}(x::Real, ::Type{T}) = reinterpret(T, _unsafe_trunc(rawtype(T), round(rawone(T)*x)))

# convert(::Type{AbstractFloat}, x::UFixed) = convert(floattype(x), x)
float(x::UFixed) = convert(floattype(x), x)

convert(::Type{BigFloat}, x::UFixed) = reinterpret(x)*(1/BigFloat(rawone(x)))
function convert{T<:AbstractFloat}(::Type{T}, x::UFixed)
    y = reinterpret(x)*(one(rawtype(x))/convert(T, rawone(x)))
    convert(T, y)  # needed for types like Float16 which promote arithmetic to Float32
end
convert(::Type{Bool}, x::UFixed) = x == zero(x) ? false : true
convert{T<:Integer}(::Type{T}, x::UFixed) = convert(T, x*(1/one(T)))
convert{Ti<:Integer}(::Type{Rational{Ti}}, x::UFixed) = convert(Ti, reinterpret(x))//convert(Ti, rawone(x))
convert(::Type{Rational}, x::UFixed) = reinterpret(x)//rawone(x)

# Traits
abs(x::UFixed) = x

(-){T<:UFixed}(x::T) = T(-reinterpret(x), 0)
(~){T<:UFixed}(x::T) = T(~reinterpret(x), 0)

+{T,f}(x::UFixed{T,f}, y::UFixed{T,f}) = UFixed{T,f}(convert(T, x.i+y.i),0)
-{T,f}(x::UFixed{T,f}, y::UFixed{T,f}) = UFixed{T,f}(convert(T, x.i-y.i),0)
*{T<:UFixed}(x::T, y::T) = convert(T,convert(floattype(T), x)*convert(floattype(T), y))
/{T<:UFixed}(x::T, y::T) = convert(T,convert(floattype(T), x)/convert(floattype(T), y))

# Comparisons
 <{T<:UFixed}(x::T, y::T) = reinterpret(x) < reinterpret(y)
<={T<:UFixed}(x::T, y::T) = reinterpret(x) <= reinterpret(y)

# Functions
trunc{T<:UFixed}(x::T) = T(div(reinterpret(x), rawone(T))*rawone(T),0)
floor{T<:UFixed}(x::T) = trunc(x)
function round{T,f}(x::UFixed{T,f})
    mask = convert(T, 1<<(f-1))
    y = trunc(x)
    return convert(T, reinterpret(x)-reinterpret(y)) & mask>0 ?
            UFixed{T,f}(y+one(UFixed{T,f})) : y
end
function ceil{T,f}(x::UFixed{T,f})
    k = 8*sizeof(T)-f
    mask = (typemax(T)<<k)>>k
    y = trunc(x)
    return convert(T, reinterpret(x)-reinterpret(y)) & (mask)>0 ?
            UFixed{T,f}(y+one(UFixed{T,f})) : y
end

trunc{T<:Integer}(::Type{T}, x::UFixed) = convert(T, div(reinterpret(x), rawone(x)))
round{T<:Integer}(::Type{T}, x::UFixed) = round(T, reinterpret(x)/rawone(x))
floor{T<:Integer}(::Type{T}, x::UFixed) = trunc(T, x)
 ceil{T<:Integer}(::Type{T}, x::UFixed) =  ceil(T, reinterpret(x)/rawone(x))

isfinite(x::UFixed) = true
isnan(x::UFixed) = false
isinf(x::UFixed) = false

bswap{f}(x::UFixed{UInt8,f}) = x
bswap(x::UFixed)  = typeof(x)(bswap(reinterpret(x)),0)

function minmax{T<:UFixed}(x::T, y::T)
    a, b = minmax(reinterpret(x), reinterpret(y))
    T(a,0), T(b,0)
end

# Iteration
# The main subtlety here is that iterating over 0x00uf8:0xffuf8 will wrap around
# unless we iterate using a wider type
@inline start{T<:UFixed}(r::StepRange{T}) = widen1(reinterpret(r.start))
@inline next{T<:UFixed}(r::StepRange{T}, i::Integer) = (T(i,0), i+reinterpret(r.step))
@inline function done{T<:UFixed}(r::StepRange{T}, i::Integer)
    i1, i2 = reinterpret(r.start), reinterpret(r.stop)
    isempty(r) | (i < min(i1, i2)) | (i > max(i1, i2))
end

function decompose(x::UFixed)
    g = gcd(reinterpret(x), rawone(x))
    div(reinterpret(x),g), 0, div(rawone(x),g)
end

# Promotions
promote_rule{T<:UFixed,Tf<:AbstractFloat}(::Type{T}, ::Type{Tf}) = promote_type(floattype(T), Tf)
promote_rule{T<:UFixed, R<:Rational}(::Type{T}, ::Type{R}) = R
function promote_rule{T<:UFixed, Ti<:Union{Signed,Unsigned}}(::Type{T}, ::Type{Ti})
    floattype(T)
end
@generated function promote_rule{T1,T2,f1,f2}(::Type{UFixed{T1,f1}}, ::Type{UFixed{T2,f2}})
    f = max(f1, f2)  # ensure we have enough precision
    T = promote_type(T1, T2)
    # make sure we have enough integer bits
    i1, i2 = 8*sizeof(T1)-f1, 8*sizeof(T2)-f2  # number of integer bits for each
    i = 8*sizeof(T)-f
    while i < max(i1, i2)
        T = widen1(T)
        i = 8*sizeof(T)-f
    end
    :(UFixed{$T,$f})
end

_unsafe_trunc{T}(::Type{T}, x::Integer) = x % T
_unsafe_trunc{T}(::Type{T}, x)          = unsafe_trunc(T, x)
