VERSION >= v"0.4.0-dev+6521" && __precompile__()

module FixedPointNumbers

using Compat

using Base: IdFun, AddFun, MulFun, reducedim_initarray

import Base: ==, <, <=, -, +, *, /, ~,
             convert, promote_rule, show, showcompact, isinteger, abs, decompose,
             isnan, isinf, isfinite,
             zero, one, typemin, typemax, realmin, realmax, eps, sizeof, reinterpret,
             trunc, round, floor, ceil, bswap,
             div, fld, rem, mod, mod1, rem1, fld1, min, max,
             start, next, done, r_promote, reducedim_init

export
    FixedPoint,
    Fixed,
    UFixed,
    Fixed16,
    UFixed8,
    UFixed10,
    UFixed12,
    UFixed14,
    UFixed16,
    # constructors
    ufixed8,
    ufixed10,
    ufixed12,
    ufixed14,
    ufixed16,
    fixed16,
    # literal constructor constants
    uf8,
    uf10,
    uf12,
    uf14,
    uf16,
    # Functions
    scaledual

# T => BaseType
# f => Number of Bytes reserved for fractional part
immutable FixedPoint{T <: Integer, f} <: Real
    i::T

    # constructor for manipulating the representation;
    # selected by passing an extra dummy argument
    FixedPoint(i::T, _) = new(i)
    FixedPoint(i::Integer,_) = new(i % T)

    FixedPoint(x) = convert(FixedPoint{T,f}, x)
end

# basic typealiases for signed and unsigned
typealias Fixed{T <: Signed, f} FixedPoint{T, f}
typealias UFixed{T <: Unsigned, f} FixedPoint{T, f}

# default provided typealiases
typealias Fixed16 Fixed{Int32, 16}
typealias UFixed8  UFixed{UInt8,8}
typealias UFixed10 UFixed{UInt16,10}
typealias UFixed12 UFixed{UInt16,12}
typealias UFixed14 UFixed{UInt16,14}
typealias UFixed16 UFixed{UInt16,16}

reinterpret(x::FixedPoint) = x.i
reinterpret{T,f}(::Type{T}, x::FixedPoint{T,f}) = x.i
reinterpret{T <: Integer,f}(::Type{FixedPoint{T,f}}, x::T) = FixedPoint{T,f}(x, 0)

rawtype{T,f}(::Type{FixedPoint{T,f}}) = T
nbitsfrac{T,f}(::Type{FixedPoint{T,f}}) = f
rawone(v) = reinterpret(one(v))

# comparisons
=={T <: FixedPoint}(x::T, y::T) = reinterpret(x) == reinterpret(y)
 <{T <: FixedPoint}(x::T, y::T) = reinterpret(x)  < reinterpret(y)
<={T <: FixedPoint}(x::T, y::T) = reinterpret(x) <= reinterpret(y)

# predicates
isinteger{T,f}(x::FixedPoint{T,f}) = (x.i&(1<<f-1)) == 0
isfinite(x::FixedPoint) = true
isnan(x::FixedPoint) = false
isinf(x::FixedPoint) = false

# traits
typemax{T<: FixedPoint}(::Type{T}) = T(typemax(rawtype(T)), 0)
typemin{T<: FixedPoint}(::Type{T}) = T(typemin(rawtype(T)), 0)
realmin{T<: FixedPoint}(::Type{T}) = typemin(T)
realmax{T<: FixedPoint}(::Type{T}) = typemax(T)
eps{T<:FixedPoint}(::Type{T}) = T(one(rawtype(T)),0)
eps{T<:FixedPoint}(::T) = eps(T)
sizeof{T<:FixedPoint}(::Type{T}) = sizeof(rawtype(T))

zero{T <: FixedPoint}(::Type{T}) = T(zero(rawtype(T)),0)
 one{T <: FixedPoint}(::Type{T}) = T(2^nbitsfrac(T)-1,0)
zero(x::FixedPoint) = zero(typeof(x))
 one(x::FixedPoint) =  one(typeof(x))


# basic operators & arithmetics
(-){T<:FixedPoint}(x::T) = T(-reinterpret(x), 0)
(~){T<:FixedPoint}(x::T) = T(~reinterpret(x), 0)
abs{T<:FixedPoint}(x::T) = T(abs(reinterpret(x)),0)

+{T<:FixedPoint}(x::T, y::T) = T(reinterpret(x) + reinterpret(y),0)
-{T<:FixedPoint}(x::T, y::T) = T(reinterpret(x) - reinterpret(y),0)

# with truncation:
# *{T<:FixedPoint}(x::T, y::T) =
#         T(Base.widemul(reinterpret(x),reinterpret(y)) >> nbitsfrac(T),0)
# with rounding up:
function *{T<:FixedPoint}(x::T, y::T)
    f = nbitsfrac(T)
    i = Base.widemul(reinterpret(x),reinterpret(y))
    T((i + convert(widen(rawtype(T)), 1) << (f-1) )>>f,0)
end

function /{T<:FixedPoint}(x::T, y::T)
    f = nbitsfrac(T)
    T(div(convert(widen(rawtype(T)), x.i) << f, y.i), 0)
end

# Conversions to FixedPoint
convert{T,f}(::Type{FixedPoint{T,f}}, x::Integer) =
    FixedPoint{T,f}(convert(T,x)<<f,0)
convert{T,f}(::Type{FixedPoint{T,f}}, x::AbstractFloat) =
    FixedPoint{T,f}(trunc(T,x) << f + round(T, rem(x,1)*(1<<f)),0)
convert{T,f}(::Type{Fixed{T,f}}, x::Rational) =
    FixedPoint{T,f}(x.num)/FixedPoint{T,f}(x.den)

# Conversions from FixedPoint
convert{T,f}(::Type{BigFloat}, x::FixedPoint{T,f}) =
    convert(BigFloat,x.i>>f) + convert(BigFloat,x.i&(1<<f - 1))/convert(BigFloat,1<<f)
# convert(::Type{BigFloat}, x::UFixed) = reinterpret(x)*(1/BigFloat(rawone(x)))

convert{TF<:AbstractFloat,T,f}(::Type{TF}, x::FixedPoint{T,f}) =
    convert(TF,x.i>>f) + convert(TF,x.i&(1<<f - 1))/convert(TF,1<<f)
# convert{T<:AbstractFloat}(::Type{T}, x::FixedPoint) = reinterpret(x)*(1/convert(T, rawone(x)))

convert{T,f}(::Type{Bool}, x::FixedPoint{T,f}) = reinterpret(x) != zero(x)

function convert{TI<:Integer, T,f}(::Type{TI}, x::FixedPoint{T,f})
    isinteger(x) || throw(InexactError())
    convert(TI, x.i>>f)
end
# convert{T<:Integer}(::Type{T}, x::UFixed) = convert(T, x*(1/one(T)))

convert{TR<:Rational,T,f}(::Type{TR}, x::FixedPoint{T,f}) =
    convert(TR, x.i>>f + (x.i&(1<<f-1))//(1<<f))
# convert{Ti<:Integer}(::Type{Rational{Ti}}, x::UFixed) = convert(Ti, reinterpret(x))//convert(Ti, rawone(x))
# convert(::Type{Rational}, x::UFixed) = reinterpret(x)//rawone(x)

# Special conversions for constructors
convert{T<:FixedPoint}(::Type{T}, x::T) = x
convert{T1<:FixedPoint}(::Type{T1}, x::FixedPoint) = reinterpret(T1, round(rawtype(T1), (rawone(T1)/rawone(x))*reinterpret(x)))
convert(::Type{UFixed16}, x::UFixed8) = reinterpret(UFixed16, convert(UInt16, 0x0101*reinterpret(x)))
convert{T<:FixedPoint}(::Type{T}, x::Real) = T(round(rawtype(T), rawone(T)*x),0)

# Constructors
ufixed8(x)  = convert(UFixed8, x)
ufixed10(x) = convert(UFixed10, x)
ufixed12(x) = convert(UFixed12, x)
ufixed14(x) = convert(UFixed14, x)
ufixed16(x) = convert(UFixed16, x)
fixed16(x) = convert(Fixed16, x)

@vectorize_1arg Real ufixed8
@vectorize_1arg Real ufixed10
@vectorize_1arg Real ufixed12
@vectorize_1arg Real ufixed14
@vectorize_1arg Real ufixed16
@vectorize_1arg Real fixed16

# Promote rules
promote_rule{T <: FixedPoint,TI<:Integer}(::Type{T}, ::Type{TI}) = T
promote_rule{T <: FixedPoint,TF<:AbstractFloat}(::Type{T}, ::Type{TF}) = TF
promote_rule{T <: FixedPoint,TR <: Rational}(::Type{T}, ::Type{TR}) = TR

# for T in UF
#     for Ti in (Int8, UInt8, Int16, UInt16, Int32, UInt32, Int64, UInt64)
#         Tp = eps(convert(Float32, typemax(Ti))) > eps(T) ? Float64 : Float32
#         @eval begin
#             promote_rule(::Type{$T}, ::Type{$Ti}) = $Tp
#         end
#     end
# end

# Math functions
trunc{T<:FixedPoint}(x::T) = T(div(reinterpret(x), rawone(T))*rawone(T),0)
floor{T<:FixedPoint}(x::T) = trunc(x)

trunc{T<:Integer}(::Type{T}, x::FixedPoint) = convert(T, div(reinterpret(x), rawone(x)))
round{T<:Integer}(::Type{T}, x::FixedPoint) = round(T, reinterpret(x)/rawone(x))
floor{T<:Integer}(::Type{T}, x::FixedPoint) = trunc(T, x)
 ceil{T<:Integer}(::Type{T}, x::FixedPoint) =  ceil(T, reinterpret(x)/rawone(x))
trunc(x::FixedPoint) = trunc(Int, x) # Already defined?
round(x::FixedPoint) = round(Int, x)
floor(x::FixedPoint) = floor(Int, x) # Already defined?
 ceil(x::FixedPoint) =  ceil(Int, x)

# for T in UF
#  f = nbitsfrac(T)
#  R = rawtype(T)
#  roundmask = convert(R, 1<<(f-1))
#  k = 8*sizeof(R)-f
#  ceilmask  = (typemax(R)<<k)>>k
#  @eval begin
#      round(x::$T) = (y = trunc(x); return convert(rawtype($T), reinterpret(x)-reinterpret(y))&$roundmask>0 ? $T(y+one($T)) : y)
#       ceil(x::$T) = (y = trunc(x); return convert(rawtype($T), reinterpret(x)-reinterpret(y))&$ceilmask >0 ? $T(y+one($T)) : y)
#  end
# end

for f in (:div, :fld, :rem, :mod, :mod1, :rem1, :fld1, :min, :max)
 @eval begin
     $f{T<:FixedPoint}(x::T, y::T) = T($f(reinterpret(x),reinterpret(y)),0)
 end
end

function minmax{T<:UFixed}(x::T, y::T)
 a, b = minmax(reinterpret(x), reinterpret(y))
 T(a,0), T(b,0)
end

# Special function
decompose{T,f}(x::FixedPoint{T,f}) = reinterpret(x), -f, 1
# function decompose(x::UFixed)
#     g = gcd(reinterpret(x), rawone(x))
#     div(reinterpret(x),g), 0, div(rawone(x),g)
# end

bswap{T <: Union{UInt8, Int8}, f}(x::FixedPoint{T,f}) = x
bswap{T <: FixedPoint}(x::T)  = T(bswap(reinterpret(x)),0)

# Promotions for reductions
const Treduce = Float64
for F in (AddFun, MulFun)
    @eval r_promote{T}(::$F, x::FixedPoint{T}) = Treduce(x)
end

reducedim_init{T<:FixedPoint}(f::IdFun, op::AddFun,
                              A::AbstractArray{T}, region) =
    reducedim_initarray(A, region, zero(Treduce))
reducedim_init{T<:FixedPoint}(f::IdFun, op::MulFun,
                              A::AbstractArray{T}, region) =
    reducedim_initarray(A, region, one(Treduce))

# When multiplying by a float, reduce two multiplies to one.
# Particularly useful for arrays.
scaledual(Tdual::Type, x) = one(Tdual), x
scaledual{Tdual<:Number}(b::Tdual, x) = b, x
@compat scaledual{T<:FixedPoint}(Tdual::Type, x::Union{T,AbstractArray{T}}) =
    convert(Tdual, 1/one(T)), reinterpret(rawtype(T), x)
@compat scaledual{Tdual<:Number, T<:FixedPoint}(b::Tdual, x::Union{T,AbstractArray{T}}) =
    convert(Tdual, b/one(T)), reinterpret(rawtype(T), x)

# Show
function show(io::IO, x::FixedPoint)
    print(io, typeof(x))
    print(io, "(")
    showcompact(io, x)
    print(io, ")")
end

const _log2_10 = 3.321928094887362
showcompact{T,f}(io::IO, x::FixedPoint{T,f}) = show(io, round(convert(Float64,x), ceil(Int,f/_log2_10)))

# Iteration
# The main subtlety here is that iterating over 0x00uf8:0xffuf8 will wrap around
# unless we iterate using a wider type
start{T<:FixedPoint}(r::StepRange{T}) = convert(typeof(reinterpret(r.start)+reinterpret(r.step)), reinterpret(r.start))
next{T<:FixedPoint}(r::StepRange{T}, i::Integer) = (T(i,0), i+reinterpret(r.step))
done{T<:FixedPoint}(r::StepRange{T}, i::Integer) = isempty(r) || (i > reinterpret(r.stop))

immutable UFixedConstructor{T,f} end
*{T,f}(n::Integer, ::UFixedConstructor{T,f}) = UFixed{T,f}(n,0)
const uf8  = UFixedConstructor{UInt8,8}()
const uf10 = UFixedConstructor{UInt16,10}()
const uf12 = UFixedConstructor{UInt16,12}()
const uf14 = UFixedConstructor{UInt16,14}()
const uf16 = UFixedConstructor{UInt16,16}()

include("deprecations.jl")

end # module
