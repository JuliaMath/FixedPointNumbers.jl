VERSION < v"0.7.0-beta2.199" && __precompile__()

module FixedPointNumbers

using Base: reducedim_initarray

import Base: ==, <, <=, -, +, *, /, ~, isapprox,
             convert, promote_rule, show, isinteger, abs, decompose,
             isnan, isinf, isfinite,
             zero, oneunit, one, typemin, typemax, floatmin, floatmax, eps, sizeof, reinterpret,
             float, trunc, round, floor, ceil, bswap,
             div, fld, rem, mod, mod1, fld1, min, max, minmax,
             start, next, done, rand
if isdefined(Base, :rem1)
    import Base: rem1
end
using Base: @pure

# T => BaseType
# f => Number of Bytes reserved for fractional part
abstract type FixedPoint{T <: Integer, f} <: Real end


export
    FixedPoint,
    Fixed,
    Normed,
# "special" typealiases
    # Q and U typealiases are exported in separate source files
# Functions
    scaledual

reinterpret(x::FixedPoint) = x.i
reinterpret(::Type{T}, x::FixedPoint{T,f}) where {T,f} = x.i

# construction using the (approximate) intended value, i.e., N0f8
*(x::Real, ::Type{X}) where {X<:FixedPoint} = X(x)

# comparison
==(x::T, y::T) where {T <: FixedPoint} = x.i == y.i
 <(x::T, y::T) where {T <: FixedPoint} = x.i  < y.i
<=(x::T, y::T) where {T <: FixedPoint} = x.i <= y.i
"""
    isapprox(x::FixedPoint, y::FixedPoint; rtol=0, atol=max(eps(x), eps(y)))

For FixedPoint numbers, the default criterion is that `x` and `y` differ by no more than `eps`, the separation between adjacent fixed-point numbers.
"""
function isapprox(x::T, y::T; rtol=0, atol=max(eps(x), eps(y))) where {T <: FixedPoint}
    maxdiff = T(atol+rtol*max(abs(x), abs(y)))
    rx, ry, rd = reinterpret(x), reinterpret(y), reinterpret(maxdiff)
    abs(signed(widen1(rx))-signed(widen1(ry))) <= rd
end
function isapprox(x::FixedPoint, y::FixedPoint; rtol=0, atol=max(eps(x), eps(y)))
    isapprox(promote(x, y)...; rtol=rtol, atol=atol)
end

# predicates
isinteger(x::FixedPoint{T,f}) where {T,f} = (x.i&(1<<f-1)) == 0

# traits
typemax(::Type{T}) where {T <: FixedPoint} = T(typemax(rawtype(T)), 0)
typemin(::Type{T}) where {T <: FixedPoint} = T(typemin(rawtype(T)), 0)
floatmin(::Type{T}) where {T <: FixedPoint} = eps(T)
floatmax(::Type{T}) where {T <: FixedPoint} = typemax(T)

widen1(::Type{Int8})   = Int16
widen1(::Type{UInt8})  = UInt16
widen1(::Type{Int16})  = Int32
widen1(::Type{UInt16}) = UInt32
widen1(::Type{Int32})  = Int64
widen1(::Type{UInt32}) = UInt64
widen1(::Type{Int64})  = Int128
widen1(::Type{UInt64}) = UInt128
widen1(::Type{UInt128}) = UInt128
widen1(x::Integer) = x % widen1(typeof(x))

const ShortInts = Union{Int8,UInt8,Int16,UInt16}

floattype(::Type{FixedPoint{T,f}}) where {T <: ShortInts,f} = Float32
floattype(::Type{FixedPoint{T,f}}) where {T <: Integer,f} = Float64
floattype(::Type{F}) where {F <: FixedPoint} = floattype(supertype(F))
floattype(x::FixedPoint) = floattype(typeof(x))

nbitsfrac(::Type{FixedPoint{T,f}}) where {T <: Integer,f} = f
nbitsfrac(::Type{F}) where {F <: FixedPoint} = nbitsfrac(supertype(F))

rawtype(::Type{FixedPoint{T,f}}) where {T <: Integer,f} = T
rawtype(::Type{F}) where {F <: FixedPoint} = rawtype(supertype(F))
rawtype(x::FixedPoint) = rawtype(typeof(x))

# This IOBuffer is used during module definition to generate typealias names
_iotypealias = IOBuffer()

# Printing. These are used to generate type-symbols, so we need them
# before we include any files.
function showtype(io::IO, ::Type{X}) where {X <: FixedPoint}
    print(io, typechar(X))
    f = nbitsfrac(X)
    m = sizeof(X)*8-f-signbits(X)
    print(io, m, 'f', f)
    io
end
function show(io::IO, x::FixedPoint{T,f}) where {T,f}
    show(io, round(convert(Float64,x), digits=ceil(Int,f/_log2_10)))
    get(io, :compact, false) || showtype(io, typeof(x))
end
const _log2_10 = 3.321928094887362

function Base.showarg(io::IO, a::Array{T}, toplevel) where {T<:FixedPoint}
    toplevel || print(io, "::")
    print(io, "Array{")
    showtype(io, T)
    print(io, ",$(ndims(a))}")
    toplevel && print(io, " with eltype ", T)
end

include("fixed.jl")
include("normed.jl")
include("deprecations.jl")
const UF = (N0f8, N6f10, N4f12, N2f14, N0f16)

eps(::Type{T}) where {T <: FixedPoint} = T(oneunit(rawtype(T)),0)
eps(::T) where {T <: FixedPoint} = eps(T)
sizeof(::Type{T}) where {T <: FixedPoint} = sizeof(rawtype(T))

# Promotions for reductions
const Treduce = Float64
if isdefined(Base, :r_promote)
    # Julia v0.6
    Base.r_promote(::typeof(+), x::FixedPoint{T}) where {T} = Treduce(x)
    Base.r_promote(::typeof(*), x::FixedPoint{T}) where {T} = Treduce(x)
    Base.reducedim_init(f::typeof(identity),
                   op::typeof(+),
                   A::AbstractArray{T}, region) where {T <: FixedPoint} =
                       Base.reducedim_initarray(A, region, zero(Treduce))
    Base.reducedim_init(f::typeof(identity),
                   op::typeof(*),
                   A::AbstractArray{T}, region) where {T <: FixedPoint} =
                       Base.reducedim_initarray(A, region, oneunit(Treduce))
else
    # Julia v0.7
    Base.add_sum(x::FixedPoint, y::FixedPoint) = Treduce(x) + Treduce(y)
    Base.reduce_empty(::typeof(Base.add_sum), ::Type{F}) where {F<:FixedPoint}  = zero(Treduce)
    Base.reduce_first(::typeof(Base.add_sum), x::FixedPoint)   = Treduce(x)
    Base.mul_prod(x::FixedPoint, y::FixedPoint) = Treduce(x) * Treduce(y)
    Base.reduce_empty(::typeof(Base.mul_prod), ::Type{F}) where {F<:FixedPoint} = one(Treduce)
    Base.reduce_first(::typeof(Base.mul_prod), x::FixedPoint)  = Treduce(x)
end


for f in (:div, :fld, :fld1)
    @eval begin
        $f(x::T, y::T) where {T <: FixedPoint} = $f(reinterpret(x),reinterpret(y))
    end
end
for f in (:rem, :mod, :mod1, :rem1, :min, :max)
    if f === :rem1 && !isdefined(Base, :rem1)
        continue
    end
    @eval begin
        $f(x::T, y::T) where {T <: FixedPoint} = T($f(reinterpret(x),reinterpret(y)),0)
    end
end

# When multiplying by a float, reduce two multiplies to one.
# Particularly useful for arrays.
scaledual(Tdual::Type, x) = oneunit(Tdual), x
scaledual(b::Tdual, x) where {Tdual <: Number} = b, x
scaledual(Tdual::Type, x::Union{T,AbstractArray{T}}) where {T <: FixedPoint} =
    convert(Tdual, 1/oneunit(T)), reinterpret(rawtype(T), x)
scaledual(b::Tdual, x::Union{T,AbstractArray{T}}) where {Tdual <: Number,T <: FixedPoint} =
    convert(Tdual, b/oneunit(T)), reinterpret(rawtype(T), x)

@noinline function throw_converterror(::Type{T}, x) where {T <: FixedPoint}
    n = 2^(8*sizeof(T))
    bitstring = sizeof(T) == 1 ? "an 8-bit" : "a $(8*sizeof(T))-bit"
    io = IOBuffer()
    show(IOContext(io, :compact=>true), typemin(T)); Tmin = String(take!(io))
    show(IOContext(io, :compact=>true), typemax(T)); Tmax = String(take!(io))
    throw(ArgumentError("$T is $bitstring type representing $n values from $Tmin to $Tmax; cannot represent $x"))
end

rand(::Type{T}) where {T <: FixedPoint} = reinterpret(T, rand(rawtype(T)))
rand(::Type{T}, sz::Dims) where {T <: FixedPoint} = reinterpret(T, rand(rawtype(T), sz))

end # module
