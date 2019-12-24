"""
    Fixed{T <: Signed, f} <: FixedPoint{T, f}

`Fixed{T,f}` maps `Signed` integers from `-2^f` to `2^f` to the range
[-1.0, 1.0]. For example, `Fixed{Int8,7}` maps `-128` to `-1.0` and `127` to
`127/128 ≈ 0.992`.

There are the typealiases for `Fixed` in the `QXfY` notation, where `Y` is
the number of fractional bits (i.e. `f`), and `X+Y+1` equals the number of
underlying bits used (`+1` means the sign bit). For example, `Q0f7` is aliased
to `Fixed{Int8,7}` and `Q3f12` is aliased to `Fixed{Int16,12}`.
"""
struct Fixed{T <: Signed, f} <: FixedPoint{T, f}
    i::T

    # constructor for manipulating the representation;
    # selected by passing an extra dummy argument
    function Fixed{T, f}(i::Integer, _) where {T, f}
        if f == bitwidth(T)
            Base.depwarn("`Fixed` reserves one bit for the sign. Support for `f=$f` with raw type `T=$T` will be removed in a future release.", :Fixed)
        end
        0 <= f <= bitwidth(T) || throw(DomainError(f, "f must be between 0 and $(bitwidth(T)-1) (i.e. the number of non-sign bits of `T=$T`)")) # TODO: change the upper limit
        new{T, f}(i % T)
    end
end

Fixed{T, f}(x::AbstractChar) where {T,f} = throw(ArgumentError("Fixed cannot be constructed from a Char"))
Fixed{T, f}(x::Complex) where {T,f} = Fixed{T, f}(convert(real(typeof(x)), x))
Fixed{T, f}(x::Base.TwicePrecision) where {T,f} = Fixed{T, f}(convert(Float64, x))
Fixed{T,f}(x::Integer) where {T,f} = Fixed{T,f}(round(T, convert(widen1(T),x)<<f),0)
Fixed{T,f}(x::AbstractFloat) where {T,f} = Fixed{T,f}(round(T, trunc(widen1(T),x)<<f + rem(x,1)*(one(widen1(T))<<f)),0)
Fixed{T,f}(x::Rational) where {T,f} = Fixed{T,f}(x.num)/Fixed{T,f}(x.den)

typechar(::Type{X}) where {X <: Fixed} = 'Q'
signbits(::Type{X}) where {X <: Fixed} = 1

for T in (Int8, Int16, Int32, Int64)
    io = IOBuffer()
    for f in 0:bitwidth(T)-1
        sym = Symbol(String(take!(showtype(io, Fixed{T,f}))))
        @eval begin
            const $sym = Fixed{$T,$f}
            export $sym
        end
    end
end

function rawone(::Type{Fixed{T,f}}) where {T, f}
    f >= bitwidth(T)-1 && throw_converterror(Fixed{T,f}, 1)
    oneunit(T) << f
end

intmask(::Fixed{T,f}) where {T, f} = -oneunit(T) << f # Signed
fracmask(x::Fixed{T,f}) where {T, f} = ~intmask(x) # Signed

# unchecked arithmetic

# with truncation:
#*(x::Fixed{T,f}, y::Fixed{T,f}) = Fixed{T,f}(Base.widemul(x.i,y.i)>>f,0)
# with rounding up:
*(x::Fixed{T,f}, y::Fixed{T,f}) where {T,f} = Fixed{T,f}((Base.widemul(x.i,y.i) + (one(widen(T)) << (f-1)))>>f,0)

/(x::Fixed{T,f}, y::Fixed{T,f}) where {T,f} = Fixed{T,f}(div(convert(widen(T), x.i) << f, y.i), 0)


# # conversions and promotions
function Fixed{T,f}(x::Fixed{T2,f2}) where {T <: Integer,T2 <: Integer,f,f2}
#    reinterpret(Fixed{T,f},T(reinterpret(x)<<(f-f2)))
    U = Fixed{T,f}
    y = round(((1<<f)/(1<<f2))*reinterpret(x))
    (typemin(T) <= y) & (y <= typemax(T)) || throw_converterror(U, x)
    reinterpret(U, _unsafe_trunc(T, y))
end

rem(x::Integer, ::Type{Fixed{T,f}}) where {T,f} = Fixed{T,f}(rem(x,T)<<f,0)
rem(x::Real,    ::Type{Fixed{T,f}}) where {T,f} = Fixed{T,f}(rem(Integer(trunc(x)),T)<<f + rem(Integer(round(rem(x,1)*(one(widen1(T))<<f))),T),0)


Base.BigFloat(x::Fixed{T,f}) where {T,f} =
    BigFloat(x.i>>f) + BigFloat(x.i&(one(widen1(T))<<f - 1))/BigFloat(one(widen1(T))<<f)
(::Type{TF})(x::Fixed{T,f}) where {TF <: AbstractFloat,T,f} =
    TF(x.i>>f) + TF(x.i&(one(widen1(T))<<f - 1))/TF(one(widen1(T))<<f)

Base.Bool(x::Fixed{T,f}) where {T,f} = x.i!=0
function Base.Integer(x::Fixed{T,f}) where {T,f}
    isinteger(x) || throw(InexactError())
    Integer(x.i>>f)
end
function (::Type{TI})(x::Fixed{T,f}) where {TI <: Integer,T,f}
    isinteger(x) || throw(InexactError())
    TI(x.i>>f)
end

(::Type{TR})(x::Fixed{T,f}) where {TR <: Rational,T,f} =
    TR(x.i>>f + (x.i&(1<<f-1))//(one(widen1(T))<<f))

function trunc(x::Fixed{T,f}) where {T, f}
    f == 0 && return x
    f == bitwidth(T) && return zero(x) # TODO: remove this line
    f == bitwidth(T) - 1 && return x.i == typemin(T) ? x : zero(x)
    t = x.i & intmask(x)
    r = x.i & fracmask(x)
    _rawone = oneunit(T) << f
    reinterpret(Fixed{T,f}, (x.i < 0) & (r != 0) ? t + _rawone : t)
end
function floor(x::Fixed{T,f}) where {T, f}
    f == bitwidth(T) && x.i < 0 && throw_converterror(Fixed{T,f}, -1) # TODO: remove this line
    Fixed{T,f}(x.i & intmask(x), 0)
end
function ceil(x::Fixed{T,f}) where {T, f}
    f == 0 && return x
    upper = typemax(T) & intmask(x)
    x.i > upper && throw_converterror(Fixed{T,f}, ceil(float(x)))
    reinterpret(Fixed{T,f}, (x.i + fracmask(x)) & intmask(x))
end
function round(x::Fixed{T,f}) where {T, f}
    f == 0 && return x
    f == bitwidth(T) && return zero(x) # TODO: remove this line
    upper = intmask(x) >>> 0x1
    lower = intmask(x) >> 0x1
    if f == bitwidth(T) - 1
        x.i > upper && throw_converterror(Fixed{T,f}, @exp2(bitwidth(T)-f-1))
        return x.i < lower ? typemin(x) : zero(x)
    end
    x.i >= upper && throw_converterror(Fixed{T,f}, @exp2(bitwidth(T)-f-1))
    y = oneunit(T) << UInt8(f - 1) + x.i
    m = oneunit(T) << UInt8(f + 1) - oneunit(T)
    z = y & intmask(x)
    reinterpret(Fixed{T,f}, z - T(y & m == rawone(x)) << f)
end

function trunc(::Type{Ti}, x::Fixed{T,f}) where {Ti <: Integer, T, f}
    f == 0 && return convert(Ti, x.i)
    f == bitwidth(T) && return zero(Ti) # TODO: remove this line
    f == bitwidth(T) - 1 && return x.i == typemin(T) ? convert(Ti, -1) : zero(Ti)
    t = x.i >> f
    r = x.i & fracmask(x)
    convert(Ti, (x.i < 0) & (r != 0) ? t + oneunit(T) : t)
end
function floor(::Type{Ti}, x::Fixed{T,f}) where {Ti <: Integer, T, f}
    f == bitwidth(T) && return x.i < 0 ? convert(Ti, -1) : zero(Ti) # TODO: remove this line
    convert(Ti, x.i >> f)
end
function ceil(::Type{Ti}, x::Fixed{T,f}) where {Ti <: Integer, T, f}
    f == bitwidth(T) && return x.i > 0 ? oneunit(Ti) : zero(Ti) # TODO: remove this line
    y = x.i + fracmask(x)
    convert(Ti, x.i >= 0 ? y >>> f : y >> f)
end
function round(::Type{Ti}, x::Fixed{T,f}) where {Ti <: Integer, T, f}
    f == 0 && return convert(Ti, x.i)
    f == bitwidth(T) && return zero(Ti) # TODO: remove this line
    upper = intmask(x) >>> 0x1
    lower = intmask(x) >> 0x1
    if f == bitwidth(T) - 1
        x.i < lower && return convert(Ti, -1)
        return x.i > upper ? oneunit(Ti) : zero(Ti)
    end
    y = oneunit(T) << UInt8(f - 1) + x.i
    m = oneunit(T) << UInt8(f + 1) - oneunit(T)
    z = x.i >= 0 ? y >>> f : y >> f
    convert(Ti, z - Ti(y & m == rawone(x)))
end

promote_rule(ft::Type{Fixed{T,f}}, ::Type{TI}) where {T,f,TI <: Integer} = Fixed{T,f}
promote_rule(::Type{Fixed{T,f}}, ::Type{TF}) where {T,f,TF <: AbstractFloat} = TF
promote_rule(::Type{Fixed{T,f}}, ::Type{Rational{TR}}) where {T,f,TR} = Rational{TR}

@generated function promote_rule(::Type{Fixed{T1,f1}}, ::Type{Fixed{T2,f2}}) where {T1,T2,f1,f2}
    f = max(f1, f2)  # ensure we have enough precision
    T = promote_type(T1, T2)
    # make sure we have enough integer bits
    i1, i2 = bitwidth(T1)-f1, bitwidth(T2)-f2  # number of integer bits for each
    i = bitwidth(T)-f
    while i < max(i1, i2)
        Tw = widen1(T)
        T == Tw && break
        T = Tw
        i = bitwidth(T)-f
    end
    :(Fixed{$T,$f})
end

# TODO: Document and check that it still does the right thing.
decompose(x::Fixed{T,f}) where {T,f} = x.i, -f, 1
