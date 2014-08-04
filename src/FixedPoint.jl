module FixedPoint

import Base: convert, promote_rule, show, showcompact, isinteger, abs,
             zero, one, typemin, typemax, realmin, realmax, eps, sizeof, reinterpret,
             trunc, round, floor, ceil, itrunc, iround, ifloor, iceil, bswap,
             div, fld, rem, mod, mod1, rem1, fld1, min, max,
             start, next, done

abstract AbstractFixed <: Real
abstract  Fixed <: AbstractFixed
abstract Ufixed <: AbstractFixed  # unsigned variant

export
    AbstractFixed,
    Fixed,
    Ufixed,
    Fixed32,
    Ufixed8,
    Ufixed10,
    Ufixed12,
    Ufixed14,
    Ufixed16,
    # literal constructor constants
    uf8,
    uf10,
    uf12,
    uf14,
    uf16,
    # Functions
    scaledual

reinterpret(x::AbstractFixed) = x.i

include("fixed32.jl")
include("ufixed.jl")

for T in tuple(Fixed32, UF...)
    R = rawtype(T)
    @eval begin
        reinterpret(::Type{$R}, x::$T) = x.i
    end
end

# When multiplying by a float, reduce two multiplies to one.
# Particularly useful for arrays.
scaledual(Tdual::Type, x) = one(Tdual), x
scaledual{Tdual<:Number}(b::Tdual, x) = b, x
scaledual{T<:AbstractFixed}(Tdual::Type, x::Union(T, AbstractArray{T})) =
    convert(Tdual, 1/one(T)), reinterpret(rawtype(T), x)
scaledual{Tdual<:Number, T<:AbstractFixed}(b::Tdual, x::Union(T, AbstractArray{T})) =
    convert(Tdual, b/one(T)), reinterpret(rawtype(T), x)

end # module
