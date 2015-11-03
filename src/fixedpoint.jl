###
# FixedPoint implements a general purpose FixedPoint number.
# The underlying storage type is given by the parameter `T`,
# and the number of fractional bits is given by `f`
# The dot is just before the fractional part so in order to represent one,
# `f` needs to be atleast smaller by 1 in the unsigned case and smaller by two
# in the signed case, then the number of bits in `T`.
#
# In a further iteration of the design `FixedPoint` should be renamed to `FixedPointPoint` and the aliase `typealias FixedPoint{T <: Signed, f} FixedPointPoint{T, f}` and `typealias UFixedPoint{T <: Unsigned, f} FixedPointPoint{T, f}` introduced.
###
immutable FixedPoint{T,f} <: AbstractFixedPoint{T,  f}
    i::T

    # constructor for manipulating the representation;
    # selected by passing an extra dummy argument
    FixedPoint(i::T, _) = new(i)
    FixedPoint(i::Integer,_) = new(i % T)

    FixedPoint(x) = convert(FixedPoint{T,f}, x)
end

typealias Fixed{T <: Signed, f} FixedPoint{T, f}
typealias Fixed16 Fixed{Int32, 16}

  rawtype{T,f}(::Type{FixedPoint{T,f}}) = T
nbitsfrac{T,f}(::Type{FixedPoint{T,f}}) = f

reinterpret{T <: Integer,f}(::Type{FixedPoint{T,f}}, x::T) = FixedPoint{T,f}(x, 0)

## predicates
isinteger{T,f}(x::FixedPoint{T,f}) = (x[]&(1<<f-1)) == 0

function one{T, f}(::Type{FixedPoint{T, f}})
    if T <: Unsigned && sizeof(T) * 8 > f ||
       T <: Signed   && sizeof(T) * 8 > (f+1)
        return FixedPoint{T, f}(2^f-1,0)
    else
        throw(DomainError())
    end
end

## basic operators & arithmetics

# with truncation:
# *{T<:FixedPoint}(x::T, y::T) =
#         T(Base.widemul(x[],y[]) >> nbitsfrac(T),0)
# with rounding up:
function *{T<:FixedPoint}(x::T, y::T)
    f = nbitsfrac(T)
    i = Base.widemul(x[],y[])
    T((i + convert(widen(rawtype(T)), 1) << (f-1) )>>f,0)
end

function /{T<:FixedPoint}(x::T, y::T)
    f = nbitsfrac(T)
    T(div(convert(widen(rawtype(T)), x[]) << f, y.i), 0)
end

## rounding
# Round towards negative infinity
trunc{T<:FixedPoint}(x::T) = T(x[] & ~(1 << nbitsfrac(T) - 1), 0)
# Round towards negative infinity
floor{T<:FixedPoint}(x::T) = trunc(x)
# Round towards positive infinity
ceil{T<:FixedPoint}(x::T) = trunc(T(x[] + 1 << (nbitsfrac(T)-1), 0))
# Round towards even
function round{T<:FixedPoint}(x::T)
    even = x[] & (1 << nbitsfrac(T)) == 0
    if even
        return floor(x)
    else
        return ceil(x)
    end
end
# Round towards positive infinity
ceil{T<:FixedPoint}(x::T) = trunc(T(x[] + 1 << (nbitsfrac(T)-1), 0))
# Round towards even
function round{T<:FixedPoint}(x::T)
    even = x[] & (1 << nbitsfrac(T)) == 0
    if even
        return floor(x)
    else
        return ceil(x)
    end
end

trunc{TI<:Integer, T <: FixedPoint}(::Type{TI}, x::T) =
    convert(TI, x[] >> nbitsfrac(T))
floor{T<:Integer}(::Type{T}, x::FixedPoint) = trunc(T, x)
ceil{T<:Integer}(::Type{T}, x::FixedPoint) = trunc(T, ceil(x))
round{T<:Integer}(::Type{T}, x::FixedPoint) = trunc(T, round(x))


## conversions and promotions
convert{T,f}(::Type{FixedPoint{T,f}}, x::Integer) = FixedPoint{T,f}(convert(T,x)<<f,0)
convert{T,f}(::Type{FixedPoint{T,f}}, x::AbstractFloat) = FixedPoint{T,f}(trunc(T,x)<<f + round(T, rem(x,1)*(1<<f)),0)
convert{T,f}(::Type{FixedPoint{T,f}}, x::Rational) = FixedPoint{T,f}(x.num)/FixedPoint{T,f}(x.den)

convert{T,f}(::Type{BigFloat}, x::FixedPoint{T,f}) =
    convert(BigFloat,x[]>>f) + convert(BigFloat,x[]&(1<<f - 1))/convert(BigFloat,1<<f)
convert{TF<:AbstractFloat,T,f}(::Type{TF}, x::FixedPoint{T,f}) =
    convert(TF,x[]>>f) + convert(TF,x[]&(1<<f - 1))/convert(TF,1<<f)

convert{T,f}(::Type{Bool}, x::FixedPoint{T,f}) = x[]!=0
function convert{TI<:Integer, T,f}(::Type{TI}, x::FixedPoint{T,f})
    isinteger(x) || throw(InexactError())
    convert(TI, x[]>>f)
end

convert{TR<:Rational,T,f}(::Type{TR}, x::FixedPoint{T,f}) =
    convert(TR, x[]>>f + (x[]&(1<<f-1))//(1<<f))

promote_rule{T,f,TI<:Integer}(ft::Type{FixedPoint{T,f}}, ::Type{TI}) = FixedPoint{T,f}
promote_rule{T,f,TF<:AbstractFloat}(::Type{FixedPoint{T,f}}, ::Type{TF}) = TF
promote_rule{T,f,TR}(::Type{FixedPoint{T,f}}, ::Type{Rational{TR}}) = Rational{TR}

decompose{T,f}(x::FixedPoint{T,f}) = x[], -f, 1

# printing
function show(io::IO, x::FixedPoint)
    print(io, typeof(x))
    print(io, "(")
    showcompact(io, x)
    print(io, ")")
end
const _log2_10 = 3.321928094887362
showcompact{T,f}(io::IO, x::FixedPoint{T,f}) = show(io, round(convert(Float64,x), ceil(Int,f/_log2_10)))
