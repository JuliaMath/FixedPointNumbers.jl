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
abstract AbstractFixedPoint{T <: Integer, f} <: Real

export
    AbstractFixedPoint,
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

getindex(x::AbstractFixedPoint) = x.i
reinterpret{T,f}(::Type{T}, x::AbstractFixedPoint{T,f}) = x[]
# reinterpret{T <: Integer,f}(::Type{AbstractFixedPoint{T,f}}, x::T) = AbstractFixedPoint{T,f}(x, 0)


# comparison
=={T <: AbstractFixedPoint}(x::T, y::T) = x[] == y[]
 <{T <: AbstractFixedPoint}(x::T, y::T) = x[]  < y[]
<={T <: AbstractFixedPoint}(x::T, y::T) = x[] <= y[]

#predicates
isfinite(x::AbstractFixedPoint) = true
isnan(x::AbstractFixedPoint) = false
isinf(x::AbstractFixedPoint) = false

#traits
typemax{T<: AbstractFixedPoint}(::Type{T}) = T(typemax(rawtype(T)), 0)
typemin{T<: AbstractFixedPoint}(::Type{T}) = T(typemin(rawtype(T)), 0)
realmin{T<: AbstractFixedPoint}(::Type{T}) = typemin(T)
realmax{T<: AbstractFixedPoint}(::Type{T}) = typemax(T)
eps{T<:AbstractFixedPoint}(::Type{T}) = T(one(rawtype(T)),0)
eps{T<:AbstractFixedPoint}(::T) = eps(T)
sizeof{T<:AbstractFixedPoint}(::Type{T}) = sizeof(rawtype(T))

zero{T <: AbstractFixedPoint}(::Type{T}) = T(zero(rawtype(T)),0)
zero{T <: AbstractFixedPoint}(x::T) = zero(T)
 one{T <: AbstractFixedPoint}(x::T) =  one(T)

# Basic operators & arithmetics
(-){T<:AbstractFixedPoint}(x::T) = T(-x[], 0)
(~){T<:AbstractFixedPoint}(x::T) = T(~x[], 0)
abs{T<:AbstractFixedPoint}(x::T) = T(abs(x[]),0)

+{T<:AbstractFixedPoint}(x::T, y::T) = T(x[] + y[],0)
-{T<:AbstractFixedPoint}(x::T, y::T) = T(x[] - y[],0)

for f in (:div, :fld, :rem, :mod, :mod1, :rem1, :fld1, :min, :max)
 @eval begin
     $f{T<:AbstractFixedPoint}(x::T, y::T) = T($f(x[],y[]),0)
 end
end

function minmax{T<:AbstractFixedPoint}(x::T, y::T)
 a, b = minmax(x[], y[])
 T(a,0), T(b,0)
end

bswap{T <: Union{UInt8, Int8}, f}(x::AbstractFixedPoint{T,f}) = x
bswap{T <: AbstractFixedPoint}(x::T)  = T(bswap(x[]),0)

include("fixedpoint.jl")
include("ufixed.jl")
include("deprecations.jl")

# Promotions for reductions
const Treduce = Float64
for F in (AddFun, MulFun)
    @eval r_promote{T}(::$F, x::AbstractFixedPoint{T}) = Treduce(x)
end

reducedim_init{T<:AbstractFixedPoint}(f::IdFun, op::AddFun,
                              A::AbstractArray{T}, region) =
    reducedim_initarray(A, region, zero(Treduce))
reducedim_init{T<:AbstractFixedPoint}(f::IdFun, op::MulFun,
                              A::AbstractArray{T}, region) =
    reducedim_initarray(A, region, one(Treduce))

# Iteration
# The main subtlety here is that iterating over 0x00uf8:0xffuf8 will wrap around
# unless we iterate using a wider type
start{T<:AbstractFixedPoint}(r::StepRange{T}) = convert(typeof(r.start[] + r.step[]), r.start[])
next{T<:AbstractFixedPoint}(r::StepRange{T}, i::Integer) = (T(i,0), i+r.step[])
done{T<:AbstractFixedPoint}(r::StepRange{T}, i::Integer) = isempty(r) || (i > r.stop[])

# When multiplying by a float, reduce two multiplies to one.
# Particularly useful for arrays.
scaledual(Tdual::Type, x) = one(Tdual), x
scaledual{Tdual<:Number}(b::Tdual, x) = b, x
scaledual{T<:AbstractFixedPoint}(Tdual::Type, x::Union{T,AbstractArray{T}}) =
    convert(Tdual, 1/one(T)), reinterpret(rawtype(T), x)
scaledual{Tdual<:Number, T<:AbstractFixedPoint}(b::Tdual, x::Union{T,AbstractArray{T}}) =
    convert(Tdual, b/one(T)), reinterpret(rawtype(T), x)

end # module
