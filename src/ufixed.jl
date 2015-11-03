# UFixed{T,f} maps UInts from 0 to 2^f-1 to the range [0.0, 1.0]
# For example, a UFixed8 maps 0x00 to 0.0 and 0xff to 1.0

immutable UFixed{T<:Unsigned,f} <: FixedPoint{T,f}
    i::T

    UFixed(i::Integer,_) = new(i%T)   # for setting by raw representation
    UFixed(x) = convert(UFixed{T,f}, x)
end

  rawtype{T,f}(::Type{UFixed{T,f}}) = T
nbitsfrac{T,f}(::Type{UFixed{T,f}}) = f

typealias UFixed8  UFixed{UInt8,8}
typealias UFixed10 UFixed{UInt16,10}
typealias UFixed12 UFixed{UInt16,12}
typealias UFixed14 UFixed{UInt16,14}
typealias UFixed16 UFixed{UInt16,16}

const UF = (UFixed8, UFixed10, UFixed12, UFixed14, UFixed16)

for (uf) in UF
    T = rawtype(uf)
    f = nbitsfrac(uf)
    @eval reinterpret(::Type{UFixed{$T,$f}}, x::$T) = UFixed{$T,$f}(x, 0)
end

# The next lines mimic the floating-point literal syntax "3.2f0"
immutable UFixedConstructor{T,f} end
*{T,f}(n::Integer, ::UFixedConstructor{T,f}) = UFixed{T,f}(n,0)
const uf8  = UFixedConstructor{UInt8,8}()
const uf10 = UFixedConstructor{UInt16,10}()
const uf12 = UFixedConstructor{UInt16,12}()
const uf14 = UFixedConstructor{UInt16,14}()
const uf16 = UFixedConstructor{UInt16,16}()

zero{T,f}(::Type{UFixed{T,f}}) = UFixed{T,f}(zero(T),0)
for uf in UF
    TT = rawtype(uf)
    f = nbitsfrac(uf)
    T = UFixed{TT,f}
    @eval begin
        one(::Type{$T}) = $T($(2^f-1),0)
    end
end
zero(x::UFixed) = zero(typeof(x))
 one(x::UFixed) =  one(typeof(x))
rawone(v) = reinterpret(one(v))

# Conversions
convert{T<:UFixed}(::Type{T}, x::T) = x
convert{T1<:UFixed}(::Type{T1}, x::UFixed) = reinterpret(T1, round(rawtype(T1), (rawone(T1)/rawone(x))*reinterpret(x)))
convert(::Type{UFixed16}, x::UFixed8) = reinterpret(UFixed16, convert(UInt16, 0x0101*reinterpret(x)))
convert{T<:UFixed}(::Type{T}, x::Real) = T(round(rawtype(T), rawone(T)*x),0)

ufixed8(x)  = convert(UFixed8, x)
ufixed10(x) = convert(UFixed10, x)
ufixed12(x) = convert(UFixed12, x)
ufixed14(x) = convert(UFixed14, x)
ufixed16(x) = convert(UFixed16, x)

@vectorize_1arg Real ufixed8
@vectorize_1arg Real ufixed10
@vectorize_1arg Real ufixed12
@vectorize_1arg Real ufixed14
@vectorize_1arg Real ufixed16


convert(::Type{BigFloat}, x::UFixed) = reinterpret(x)*(1/BigFloat(rawone(x)))
convert{T<:AbstractFloat}(::Type{T}, x::UFixed) = reinterpret(x)*(1/convert(T, rawone(x)))
convert(::Type{Bool}, x::UFixed) = x == zero(x) ? false : true
convert{T<:Integer}(::Type{T}, x::UFixed) = convert(T, x*(1/one(T)))
convert{Ti<:Integer}(::Type{Rational{Ti}}, x::UFixed) = convert(Ti, reinterpret(x))//convert(Ti, rawone(x))
convert(::Type{Rational}, x::UFixed) = reinterpret(x)//rawone(x)

# Traits
eps{T<:UFixed}(::Type{T}) = T(one(rawtype(T)),0)
eps{T<:UFixed}(::T) = eps(T)
sizeof{T<:UFixed}(::Type{T}) = sizeof(rawtype(T))
abs(x::UFixed) = x

# Arithmetic
(-){T<:UFixed}(x::T) = T(-reinterpret(x), 0)
(~){T<:UFixed}(x::T) = T(~reinterpret(x), 0)

+{T,f}(x::UFixed{T,f}, y::UFixed{T,f}) = UFixed{T,f}(convert(T, x.i+y.i),0)
-{T,f}(x::UFixed{T,f}, y::UFixed{T,f}) = UFixed{T,f}(convert(T, x.i-y.i),0)
*{T<:UFixed}(x::T, y::T) = convert(T,convert(Float32, x)*convert(Float32, y))
/{T<:UFixed}(x::T, y::T) = convert(T,convert(Float32, x)/convert(Float32, y))

# Comparisons
 <{T<:UFixed}(x::T, y::T) = reinterpret(x) < reinterpret(y)
<={T<:UFixed}(x::T, y::T) = reinterpret(x) < reinterpret(y)

# Functions
trunc{T<:UFixed}(x::T) = T(div(reinterpret(x), rawone(T))*rawone(T),0)
floor{T<:UFixed}(x::T) = trunc(x)
for T in UF
    f = nbitsfrac(T)
    R = rawtype(T)
    roundmask = convert(R, 1<<(f-1))
    k = 8*sizeof(R)-f
    ceilmask  = (typemax(R)<<k)>>k
    @eval begin
        round(x::$T) = (y = trunc(x); return convert(rawtype($T), reinterpret(x)-reinterpret(y))&$roundmask>0 ? $T(y+one($T)) : y)
         ceil(x::$T) = (y = trunc(x); return convert(rawtype($T), reinterpret(x)-reinterpret(y))&$ceilmask >0 ? $T(y+one($T)) : y)
    end
end

trunc{T<:Integer}(::Type{T}, x::UFixed) = convert(T, div(reinterpret(x), rawone(x)))
round{T<:Integer}(::Type{T}, x::UFixed) = round(T, reinterpret(x)/rawone(x))
floor{T<:Integer}(::Type{T}, x::UFixed) = trunc(T, x)
 ceil{T<:Integer}(::Type{T}, x::UFixed) =  ceil(T, reinterpret(x)/rawone(x))
trunc(x::UFixed) = trunc(Int, x)
round(x::UFixed) = round(Int, x)
floor(x::UFixed) = floor(Int, x)
 ceil(x::UFixed) =  ceil(Int, x)

isfinite(x::UFixed) = true
isnan(x::UFixed) = false
isinf(x::UFixed) = false

bswap{f}(x::UFixed{UInt8,f}) = x
bswap(x::UFixed)  = typeof(x)(bswap(reinterpret(x)),0)

for f in (:div, :fld, :rem, :mod, :mod1, :rem1, :fld1, :min, :max)
    @eval begin
        $f{T<:UFixed}(x::T, y::T) = T($f(reinterpret(x),reinterpret(y)),0)
    end
end
function minmax{T<:UFixed}(x::T, y::T)
    a, b = minmax(reinterpret(x), reinterpret(y))
    T(a,0), T(b,0)
end

# Iteration
# The main subtlety here is that iterating over 0x00uf8:0xffuf8 will wrap around
# unless we iterate using a wider type
if VERSION < v"0.3-"
    start{T<:UFixed}(r::Range{T}) = convert(typeof(reinterpret(r.start)+reinterpret(r.step)), reinterpret(r.start))
    next{T<:UFixed}(r::Range{T}, i::Integer) = (T(i,0), i+reinterpret(r.step))
    done{T<:UFixed}(r::Range{T}, i::Integer) = isempty(r) || (i > r.len)
else
    start{T<:UFixed}(r::StepRange{T}) = convert(typeof(reinterpret(r.start)+reinterpret(r.step)), reinterpret(r.start))
    next{T<:UFixed}(r::StepRange{T}, i::Integer) = (T(i,0), i+reinterpret(r.step))
    done{T<:UFixed}(r::StepRange{T}, i::Integer) = isempty(r) || (i > reinterpret(r.stop))
end

function decompose(x::UFixed)
    g = gcd(reinterpret(x), rawone(x))
    div(reinterpret(x),g), 0, div(rawone(x),g)
end

# Promotions
for T in UF
    @eval begin
        promote_rule(::Type{$T}, ::Type{Float32}) = Float32
        promote_rule(::Type{$T}, ::Type{Float64}) = Float64
        promote_rule{TR<:Rational}(::Type{$T}, ::Type{TR}) = TR
    end
    for Ti in (Int8, UInt8, Int16, UInt16, Int32, UInt32, Int64, UInt64)
        Tp = eps(convert(Float32, typemax(Ti))) > eps(T) ? Float64 : Float32
        @eval begin
            promote_rule(::Type{$T}, ::Type{$Ti}) = $Tp
        end
    end
end

# Show
function show{T,f}(io::IO, x::UFixed{T,f})
    print(io, "UFixed", f)
    print(io, "(")
    showcompact(io, x)
    print(io, ")")
end
showcompact{T,f}(io::IO, x::UFixed{T,f}) = show(io, round(convert(Float64,x), ceil(Int,f/_log2_10)))
