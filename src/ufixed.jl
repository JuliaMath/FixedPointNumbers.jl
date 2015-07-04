# UfixedBase{T,f} maps Uints from 0 to 2^f-1 to the range [0.0, 1.0]
# For example, a Ufixed8 maps 0x00 to 0.0 and 0xff to 1.0

immutable UfixedBase{T<:Unsigned,f} <: Ufixed
    i::T

    UfixedBase(i::Integer,_) = new(i)   # for setting by raw representation
    UfixedBase(x) = convert(UfixedBase{T,f}, x)
end

typealias Ufixed8  UfixedBase{Uint8,8}
typealias Ufixed10 UfixedBase{Uint16,10}
typealias Ufixed12 UfixedBase{Uint16,12}
typealias Ufixed14 UfixedBase{Uint16,14}
typealias Ufixed16 UfixedBase{Uint16,16}

const UF = (Ufixed8, Ufixed10, Ufixed12, Ufixed14, Ufixed16)

  rawtype{T,f}(::Type{UfixedBase{T,f}}) = T
nbitsfrac{T,f}(::Type{UfixedBase{T,f}}) = f

reinterpret(::Type{Ufixed8}, x::Uint8) = UfixedBase{Uint8,8}(x,0)
for (T,f) in ((Ufixed10,10),(Ufixed12,12),(Ufixed14,14),(Ufixed16,16))
    @eval reinterpret(::Type{$T}, x::Uint16) = UfixedBase{Uint16,$f}(x, 0)
end


# The next lines mimic the floating-point literal syntax "3.2f0"
immutable UfixedConstructor{T,f} end
*{T,f}(n::Integer, ::UfixedConstructor{T,f}) = UfixedBase{T,f}(n,0)
const uf8  = UfixedConstructor{Uint8,8}()
const uf10 = UfixedConstructor{Uint16,10}()
const uf12 = UfixedConstructor{Uint16,12}()
const uf14 = UfixedConstructor{Uint16,14}()
const uf16 = UfixedConstructor{Uint16,16}()

zero{T,f}(::Type{UfixedBase{T,f}}) = UfixedBase{T,f}(zero(T),0)
for T in UF
    f = nbitsfrac(T)
    @eval begin
        one(::Type{$T}) = UfixedBase{$(rawtype(T)),$f}($(2^f-1),0)
    end
end
zero(x::Ufixed) = zero(typeof(x))
 one(x::Ufixed) =  one(typeof(x))
rawone(v) = reinterpret(one(v))

# Conversions
convert{T<:Ufixed}(::Type{T}, x::T)    = x
convert{T1<:Ufixed}(::Type{T1}, x::Ufixed)  = reinterpret(T1, round(rawtype(T1), (rawone(T1)/rawone(x))*reinterpret(x)))
convert(::Type{Ufixed16}, x::Ufixed8)  = reinterpret(Ufixed16, convert(Uint16, 0x0101*reinterpret(x)))
convert{T<:Ufixed}(::Type{T}, x::Real) = T(round(rawtype(T), rawone(T)*x),0)

ufixed8(x)  = convert(Ufixed8, x)
ufixed10(x) = convert(Ufixed10, x)
ufixed12(x) = convert(Ufixed12, x)
ufixed14(x) = convert(Ufixed14, x)
ufixed16(x) = convert(Ufixed16, x)

@vectorize_1arg Real ufixed8
@vectorize_1arg Real ufixed10
@vectorize_1arg Real ufixed12
@vectorize_1arg Real ufixed14
@vectorize_1arg Real ufixed16


convert(::Type{BigFloat}, x::Ufixed) = reinterpret(x)*(1/BigFloat(rawone(x)))
convert{T<:FloatingPoint}(::Type{T}, x::Ufixed) = reinterpret(x)*(1/convert(T, rawone(x)))
convert(::Type{Bool}, x::Ufixed) = x == zero(x) ? false : true
convert{T<:Integer}(::Type{T}, x::Ufixed) = convert(T, x*(1/one(T)))
convert{Ti<:Integer}(::Type{Rational{Ti}}, x::Ufixed) = convert(Ti, reinterpret(x))//convert(Ti, rawone(x))
convert(::Type{Rational}, x::Ufixed) = reinterpret(x)//rawone(x)

# Traits
typemin{T<:Ufixed}(::Type{T}) = zero(T)
typemax{T<:Ufixed}(::Type{T}) = T(typemax(rawtype(T)),0)
realmin{T<:Ufixed}(::Type{T}) = typemin(T)
realmax{T<:Ufixed}(::Type{T}) = typemax(T)
eps{T<:Ufixed}(::Type{T}) = T(one(rawtype(T)),0)
eps{T<:Ufixed}(::T) = eps(T)
sizeof{T<:Ufixed}(::Type{T}) = sizeof(rawtype(T))
abs(x::Ufixed) = x

# Arithmetic
+{T,f}(x::UfixedBase{T,f}, y::UfixedBase{T,f}) = convert(Float32, x)+convert(Float32, y) # UfixedBase{T,f}(convert(T, reinterpret(x)+reinterpret(y)),0)
-{T,f}(x::UfixedBase{T,f}, y::UfixedBase{T,f}) = convert(Float32, x)-convert(Float32, y) # UfixedBase{T,f}(convert(T, reinterpret(x)-reinterpret(y)),0)
*{T,f}(x::UfixedBase{T,f}, y::UfixedBase{T,f}) = convert(Float32, x)*convert(Float32, y)
/(x::Ufixed, y::Ufixed) = convert(Float32, x)/convert(Float32, y)

# Comparisons
 <{T<:Ufixed}(x::T, y::T) = reinterpret(x) <  reinterpret(y)
<={T<:Ufixed}(x::T, y::T) = reinterpret(x) <  reinterpret(y)

# Functions
trunc{T<:Ufixed}(x::T) = T(div(reinterpret(x), rawone(T))*rawone(T),0)
floor{T<:Ufixed}(x::T) = trunc(x)
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

trunc{T<:Integer}(::Type{T}, x::Ufixed) = convert(T, div(reinterpret(x), rawone(x)))
round{T<:Integer}(::Type{T}, x::Ufixed) = round(T, reinterpret(x)/rawone(x))
floor{T<:Integer}(::Type{T}, x::Ufixed) = trunc(T, x)
 ceil{T<:Integer}(::Type{T}, x::Ufixed) =  ceil(T, reinterpret(x)/rawone(x))
trunc(x::Ufixed) = trunc(Int, x)
round(x::Ufixed) = round(Int, x)
floor(x::Ufixed) = floor(Int, x)
 ceil(x::Ufixed) =  ceil(Int, x)

isfinite(x::Ufixed) = true
isnan(x::Ufixed) = false
isinf(x::Ufixed) = false

bswap{f}(x::UfixedBase{Uint8,f}) = x
bswap(x::Ufixed)  = typeof(x)(bswap(reinterpret(x)),0)

for f in (:div, :fld, :rem, :mod, :mod1, :rem1, :fld1, :min, :max)
    @eval begin
        $f{T<:Ufixed}(x::T, y::T) = T($f(reinterpret(x),reinterpret(y)),0)
    end
end
function minmax{T<:Ufixed}(x::T, y::T)
    a, b = minmax(reinterpret(x), reinterpret(y))
    T(a,0), T(b,0)
end

# Iteration
# The main subtlety here is that iterating over 0x00uf8:0xffuf8 will wrap around
# unless we iterate using a wider type
if VERSION < v"0.3-"
    start{T<:Ufixed}(r::Range{T}) = convert(typeof(reinterpret(r.start)+reinterpret(r.step)), reinterpret(r.start))
    next{T<:Ufixed}(r::Range{T}, i::Integer) = (T(i,0), i+reinterpret(r.step))
    done{T<:Ufixed}(r::Range{T}, i::Integer) = isempty(r) || (i > r.len)
else
    start{T<:Ufixed}(r::StepRange{T}) = convert(typeof(reinterpret(r.start)+reinterpret(r.step)), reinterpret(r.start))
    next{T<:Ufixed}(r::StepRange{T}, i::Integer) = (T(i,0), i+reinterpret(r.step))
    done{T<:Ufixed}(r::StepRange{T}, i::Integer) = isempty(r) || (i > reinterpret(r.stop))
end

function decompose(x::Ufixed)
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
    for Ti in (Int8, Uint8, Int16, Uint16, Int32, Uint32, Int64, Uint64)
        Tp = eps(convert(Float32, typemax(Ti))) > eps(T) ? Float64 : Float32
        @eval begin
            promote_rule(::Type{$T}, ::Type{$Ti}) = $Tp
        end
    end
end

# Show
function show{T,f}(io::IO, x::UfixedBase{T,f})
    print(io, "Ufixed", f)
    print(io, "(")
    showcompact(io, x)
    print(io, ")")
end
showcompact{T,f}(io::IO, x::UfixedBase{T,f}) = show(io, round(convert(Float64,x), ceil(Int,f/_log2_10)))
