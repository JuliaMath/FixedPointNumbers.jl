__precompile__()

module FixedPointNumbers

using Base: IdFun, AddFun, MulFun, reducedim_initarray

import Base: ==, <, <=, -, +, *, /, ~,
             convert, promote_rule, show, showcompact, isinteger, abs, decompose,
             isnan, isinf, isfinite,
             zero, one, typemin, typemax, realmin, realmax, eps, sizeof, reinterpret, getindex,
             trunc, round, floor, ceil, bswap,
             div, fld, rem, mod, mod1, rem1, fld1, min, max,
             start, next, done, r_promote, reducedim_init
# T => BaseType
# f => Number of Bytes reserved for fractional part
abstract FixedPoint{T <: Integer, f} <: Real

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
    # literal constructor constants
    uf8,
    uf10,
    uf12,
    uf14,
    uf16,
    # Functions
    scaledual

getindex(x::FixedPoint) = x.i
reinterpret{T,f}(::Type{T}, x::FixedPoint{T,f}) = x[]
# reinterpret{T <: Integer,f}(::Type{FixedPoint{T,f}}, x::T) = FixedPoint{T,f}(x, 0)


# comparison
=={T <: FixedPoint}(x::T, y::T) = x[] == y[]
 <{T <: FixedPoint}(x::T, y::T) = x[]  < y[]
<={T <: FixedPoint}(x::T, y::T) = x[] <= y[]

#predicates
isfinite(x::FixedPoint) = true
isnan(x::FixedPoint) = false
isinf(x::FixedPoint) = false

#traits
typemax{T<: FixedPoint}(::Type{T}) = T(typemax(rawtype(T)), 0)
typemin{T<: FixedPoint}(::Type{T}) = T(typemin(rawtype(T)), 0)
realmin{T<: FixedPoint}(::Type{T}) = typemin(T)
realmax{T<: FixedPoint}(::Type{T}) = typemax(T)
eps{T<:FixedPoint}(::Type{T}) = T(one(rawtype(T)),0)
eps{T<:FixedPoint}(::T) = eps(T)
sizeof{T<:FixedPoint}(::Type{T}) = sizeof(rawtype(T))

zero{T <: FixedPoint}(::Type{T}) = T(zero(rawtype(T)),0)
zero{T <: FixedPoint}(x::T) = zero(T)
 one{T <: FixedPoint}(x::T) =  one(T)

# Basic operators & arithmetics
(-){T<:FixedPoint}(x::T) = T(-x[], 0)
(~){T<:FixedPoint}(x::T) = T(~x[], 0)
abs{T<:FixedPoint}(x::T) = T(abs(x[]),0)

+{T<:FixedPoint}(x::T, y::T) = T(x[] + y[],0)
-{T<:FixedPoint}(x::T, y::T) = T(x[] - y[],0)

for f in (:div, :fld, :rem, :mod, :mod1, :rem1, :fld1, :min, :max)
 @eval begin
     $f{T<:FixedPoint}(x::T, y::T) = T($f(x[],y[]),0)
 end
end

function minmax{T<:FixedPoint}(x::T, y::T)
 a, b = minmax(x[], y[])
 T(a,0), T(b,0)
end

bswap{T <: Union{UInt8, Int8}, f}(x::FixedPoint{T,f}) = x
bswap{T <: FixedPoint}(x::T)  = T(bswap(x[]),0)

include("fixed.jl")
include("ufixed.jl")
include("deprecations.jl")

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

# Iteration
# The main subtlety here is that iterating over 0x00uf8:0xffuf8 will wrap around
# unless we iterate using a wider type
start{T<:FixedPoint}(r::StepRange{T}) = convert(typeof(r.start[] + r.step[]), r.start[])
next{T<:FixedPoint}(r::StepRange{T}, i::Integer) = (T(i,0), i+r.step[])
done{T<:FixedPoint}(r::StepRange{T}, i::Integer) = isempty(r) || (i > r.stop[])

# When multiplying by a float, reduce two multiplies to one.
# Particularly useful for arrays.
scaledual(Tdual::Type, x) = one(Tdual), x
scaledual{Tdual<:Number}(b::Tdual, x) = b, x
scaledual{T<:FixedPoint}(Tdual::Type, x::Union{T,AbstractArray{T}}) =
    convert(Tdual, 1/one(T)), reinterpret(rawtype(T), x)
scaledual{Tdual<:Number, T<:FixedPoint}(b::Tdual, x::Union{T,AbstractArray{T}}) =
    convert(Tdual, b/one(T)), reinterpret(rawtype(T), x)

end # module
