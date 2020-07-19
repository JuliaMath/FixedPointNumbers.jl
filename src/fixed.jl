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

typechar(::Type{X}) where {X <: Fixed} = 'Q'

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

# constructor-style conversions
function _convert(::Type{F}, x::Fixed{T2,f2}) where {T, T2, f, f2, F <: Fixed{T,f}}
    y = round(@exp2(f-f2) * reinterpret(x))
    (typemin(T) <= y) & (y <= typemax(T)) || throw_converterror(F, x)
    reinterpret(F, _unsafe_trunc(T, y))
end

function _convert(::Type{F}, x::Integer) where {T, f, F <: Fixed{T,f}}
    if ((typemin(T) >> f) <= x) & (x <= (typemax(T) >> f))
        reinterpret(F, _unsafe_trunc(T, x) << f)
    else
        throw_converterror(F, x)
    end
end

function _convert(::Type{F}, x::AbstractFloat) where {T, f, F <: Fixed{T,f}}
    bigx = big(x)
    bmin = BigFloat(typemin(F)) - @exp2(-f-1)
    bmax = BigFloat(typemax(F)) + @exp2(-f-1)
    if bmin <= bigx < bmax
        reinterpret(F, round(T, bigx * @exp2(f)))
    else
        throw_converterror(F, x)
    end
end

_convert(::Type{F}, x::Float16) where {T, f, F <: Fixed{T,f}} = F(Float32(x))

function _convert(::Type{F}, x::Union{Float32, Float64}) where {T, f, F <: Fixed{T,f}}
    Tf = typeof(x)
    if Tf(typemin(F) - @exp2(-f-1)) <= x < Tf(typemax(F) + @exp2(-f-1))
        reinterpret(F, round(T, x * @exp2(f)))
    else
        throw_converterror(F, x)
    end
end

function _convert(::Type{F}, x::Rational) where {T, f, F <: Fixed{T,f}}
    xmin = widemul(denominator(x), widen1(T)(typemin(T)) << 0x1 - 0x1)
    xmax = widemul(denominator(x), oneunit(widen1(T)) << bitwidth(T) - 0x1)
    if xmin <= (widen1(numerator(x)) << UInt8(f + 1)) < xmax
        reinterpret(F, round(T, convert(floattype(T), x) * @exp2(f)))
    else
        throw_converterror(F, x)
    end
end

# unchecked arithmetic

# with truncation:
#*(x::Fixed{T,f}, y::Fixed{T,f}) = Fixed{T,f}(Base.widemul(x.i,y.i)>>f,0)
# with rounding up:
*(x::Fixed{T,f}, y::Fixed{T,f}) where {T,f} = Fixed{T,f}((Base.widemul(x.i,y.i) + (one(widen(T)) << (f-1)))>>f,0)

/(x::Fixed{T,f}, y::Fixed{T,f}) where {T,f} = Fixed{T,f}(div(convert(widen(T), x.i) << f, y.i), 0)


rem(x::Integer, ::Type{Fixed{T,f}}) where {T,f} = Fixed{T,f}(rem(x,T)<<f,0)
rem(x::Real,    ::Type{Fixed{T,f}}) where {T,f} = Fixed{T,f}(rem(Integer(trunc(x)),T)<<f + rem(Integer(round(rem(x,1)*(one(widen1(T))<<f))),T),0)


(::Type{Tf})(x::Fixed{T,f}) where {Tf <: AbstractFloat, T, f} = Tf(Tf(x.i) * Tf(@exp2(-f)))
Base.Float16(x::Fixed{T,f}) where {T, f} = Float16(Float32(x))
Base.Float32(x::Fixed{T,f}) where {T, f} = Float32(x.i) * Float32(@exp2(-f))
Base.Float64(x::Fixed{T,f}) where {T, f} = Float64(x.i) * @exp2(-f)

function Base.Rational(x::Fixed{T,f}) where {T, f}
    f < bitwidth(T)-1 ? x.i//rawone(x) : x.i//(one(widen1(T))<<f)
end

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

# Range construction
Base.unitrange_last(start::F, stop::F) where {F<:Fixed} =
    stop >= start ? convert(F, start+floor(stop-start)) : convert(F, start+F(-1))

# Range lengths
length(r::AbstractUnitRange{F}) where {F <: Fixed{<:SShorterThanInt,f}} where {f} =
    ((Int(reinterpret(last(r))) - Int(reinterpret(first(r)))) >> f) + 1
length(r::AbstractUnitRange{F}) where {F <: Fixed{T}} where {T <: Signed} =
    checked_add(checked_sub(floor(T, last(r)), floor(T, first(r))), oneunit(T))

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
