__precompile__()

module FixedPointNumbers

using Base: IdFun, AddFun, MulFun, reducedim_initarray

import Base: ==, <, <=, -, +, *, /, ~,
             convert, promote_rule, show, showcompact, isinteger, abs, decompose,
             isnan, isinf, isfinite,
             zero, one, typemin, typemax, realmin, realmax, eps, sizeof, reinterpret,
             trunc, round, floor, ceil, bswap,
             div, fld, rem, mod, mod1, rem1, fld1, min, max, minmax,
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

reinterpret(x::FixedPoint) = x.i

# comparison
=={T <: FixedPoint}(x::T, y::T) = x.i == y.i
 <{T <: FixedPoint}(x::T, y::T) = x.i  < y.i
<={T <: FixedPoint}(x::T, y::T) = x.i <= y.i

# predicates
isinteger{T,f}(x::FixedPoint{T,f}) = (x.i&(1<<f-1)) == 0

typemax{T<: FixedPoint}(::Type{T}) = T(typemax(rawtype(T)), 0)
typemin{T<: FixedPoint}(::Type{T}) = T(typemin(rawtype(T)), 0)
realmin{T<: FixedPoint}(::Type{T}) = typemin(T)
realmax{T<: FixedPoint}(::Type{T}) = typemax(T)

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

# TODO: rewrite this by @generated
for T in tuple(Fixed16, UF...)
    R = rawtype(T)
    @eval begin
        reinterpret(::Type{$R}, x::$T) = x.i
    end
end

# When multiplying by a float, reduce two multiplies to one.
# Particularly useful for arrays.
scaledual(Tdual::Type, x) = one(Tdual), x
scaledual{Tdual<:Number}(b::Tdual, x) = b, x
scaledual{T<:FixedPoint}(Tdual::Type, x::Union{T,AbstractArray{T}}) =
    convert(Tdual, 1/one(T)), reinterpret(rawtype(T), x)
scaledual{Tdual<:Number, T<:FixedPoint}(b::Tdual, x::Union{T,AbstractArray{T}}) =
    convert(Tdual, b/one(T)), reinterpret(rawtype(T), x)

end # module
