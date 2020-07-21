"""
    Normed{T <: Unsigned, f} <: FixedPoint{T, f}

`Normed{T,f}` maps `Unsigned` integers from `0` to `2^f-1` to the range
[0.0, 1.0]. For example, `Normed{UInt8,8}` maps `0x00` to `0.0` and `0xff` to
`1.0`.

There are the typealiases for `Normed` in the `NXfY` notation, where `Y` is
the number of fractional bits (i.e. `f`), and `X+Y` equals the number of
underlying bits used. For example, `N0f8` is aliased to `Normed{UInt8,8}` and
`N4f12` is aliased to `Normed{UInt16,12}`.
"""
struct Normed{T <: Unsigned, f} <: FixedPoint{T, f}
    i::T

    function Normed{T, f}(i::Integer, _) where {T, f}
        1 <= f <= bitwidth(T) || throw(DomainError(f, "f must be between 1 and $(bitwidth(T)) (i.e. the number of bits of `T=$T`)"))
        new{T, f}(i % T)
    end
end

typechar(::Type{X}) where {X <: Normed} = 'N'

for T in (UInt8, UInt16, UInt32, UInt64)
    io = IOBuffer()
    for f in 1:bitwidth(T)
        sym = Symbol(String(take!(showtype(io, Normed{T,f}))))
        @eval begin
            const $sym = Normed{$T,$f}
            export $sym
        end
    end
end

function rawone(::Type{Normed{T,f}}) where {T <: Unsigned, f}
    typemax(T) >> (bitwidth(T) - f)
end

# constructor-style conversions
function _convert(::Type{N}, x::Normed{T2,f}) where {T, T2, f, N <: Normed{T,f}}
    reinterpret(N, convert(T, x.i)) # TODO: input range checking
end

function _convert(::Type{N}, x::Normed{T2,f2}) where {T, T2, f, f2, N <: Normed{T,f}}
    y = round((rawone(N)/rawone(x))*reinterpret(x))
    (0 <= y) & (y <= typemax(T)) || throw_converterror(N, x)
    reinterpret(N, _unsafe_trunc(T, y))
end

function _convert(::Type{N}, x::Normed{UInt8,8}) where {N <: Normed{UInt16,16}} # TODO: generalization
    reinterpret(N0f16, convert(UInt16, 0x0101*reinterpret(x)))
end

function _convert(::Type{N}, x::Real) where {T, f, N <: Normed{T,f}}
    if T == UInt128 # for UInt128, we can't widen
        # the upper limit is not exact
        (0 <= x) & (x <= (typemax(T)/rawone(N))) || throw_converterror(N, x)
        y = round(rawone(N)*x)
    else
        y = round(widen1(rawone(N))*x)
        (0 <= y) & (y <= typemax(T)) || throw_converterror(N, x)
    end
    reinterpret(N, _unsafe_trunc(T, y))
end
# Prevent overflow (https://discourse.julialang.org/t/saving-greater-than-8-bit-images/6057)
function _convert(::Type{N}, x::Float16) where {T, f, N <: Normed{T,f}}
    if Float16(typemax(T)/rawone(N)) > Float32(typemax(T)/rawone(N))
        x == Float16(typemax(T)/rawone(N)) && return typemax(N)
    end
    return _convert(N, Float32(x))
end
function _convert(::Type{N}, x::Tf) where {T, f, N <: Normed{T,f}, Tf <: Union{Float32, Float64}}
    if T === UInt128 && f == 53
        0 <= x <= Tf(3.777893186295717e22) || throw_converterror(N, x)
    else
        0 <= x <= Tf((typemax(T)-rawone(N))/rawone(N)+1) || throw_converterror(N, x)
    end

    f == 1 && x == Tf(typemax(N)) && return typemax(N)
    if f <= (significand_bits(Tf) + 1) && bitwidth(T) < significand_bits(Tf)
        return reinterpret(N, unsafe_trunc(T, round(rawone(N) * x)))
    end
    # cf. the implementation of `frexp`
    Tw = f < bitwidth(T) ? T : widen1(T)
    bits = bitwidth(Tw) - 1
    xu = reinterpret(Unsigned, x)
    k = Int(xu >> significand_bits(Tf))
    k == 0 && return zero(N) # neglect subnormal numbers
    significand = xu | (oneunit(xu) << significand_bits(Tf))
    yh = unsafe_trunc(Tw, significand) << (bits - significand_bits(Tf))
    ex = exponent_bias(Tf) - k + bits - f
    yi = bits >= f ? yh - (yh >> f) : yh
    if ex <= 0
        ex == 0 && return reinterpret(N, unsafe_trunc(T, yi))
        ex != -1 || signbit(signed(yi)) && return typemax(N)
        return reinterpret(N, unsafe_trunc(T, yi + yi))
    end
    ex > bits && return reinterpret(N, ex == bits + 1 ? oneunit(T) : zero(T))
    yi += oneunit(Tw)<<((ex - 1) & bits) # RoundNearestTiesUp
    return reinterpret(N, unsafe_trunc(T, yi >> (ex & bits)))
end

function _convert(::Type{N}, x::Rational) where {T, f, N <: Normed{T,f}}
    if 0 <= x <= Rational(typemax(N))
        reinterpret(N, round(T, convert(floattype(T), x) * rawone(N)))
    else
        throw_converterror(N, x)
    end
end

rem(x::N, ::Type{N}) where {N <: Normed} = x
rem(x::Normed, ::Type{N}) where {T, N <: Normed{T}} = reinterpret(N, _unsafe_trunc(T, round((rawone(N)/rawone(x))*reinterpret(x))))
rem(x::Real, ::Type{N}) where {T, N <: Normed{T}} = reinterpret(N, _unsafe_trunc(T, round(rawone(N)*x)))
rem(x::Float16, ::Type{N}) where {N <: Normed} = rem(Float32(x), N)  # avoid overflow
# Float32 and Float64 cannot exactly represent `rawone(N)` with `f` greater than
# the number of their significand bits, resulting in rounding errors (issue #150).
# So, we use another strategy for the large `f`s explained in:
# https://github.com/JuliaMath/FixedPointNumbers.jl/pull/166#issuecomment-574135643
function rem(x::Float32, ::Type{N}) where {f, N <: Normed{UInt32,f}}
    f <= 24 && return reinterpret(N, _unsafe_trunc(UInt32, round(rawone(N) * x)))
    r = _unsafe_trunc(UInt32, round(x * @f32(0x1p24)))
    reinterpret(N, r << UInt8(f - 24) - unsigned(signed(r) >> 0x18))
end
function rem(x::Float64, ::Type{N}) where {f, N <: Normed{UInt64,f}}
    f <= 53 && return reinterpret(N, _unsafe_trunc(UInt64, round(rawone(N) * x)))
    r = _unsafe_trunc(UInt64, round(x * 0x1p53))
    reinterpret(N, r << UInt8(f - 53) - unsigned(signed(r) >> 0x35))
end


function (::Type{T})(x::Normed) where {T <: AbstractFloat}
    # The following optimization for constant division may cause rounding errors.
    # y = reinterpret(x)*(one(rawtype(x))/convert(T, rawone(x)))
    # Therefore, we use a simple form here.
    # If you prefer speed over accuracy, consider using `scaledual` instead.
    y = reinterpret(x) / convert(promote_type(T, floattype(x)), rawone(x))
    convert(T, y)  # needed for types like Float16 which promote arithmetic to Float32
end

function Base.Float16(x::Normed{Ti,f}) where {Ti <: Union{UInt8, UInt16, UInt32}, f}
    f == 1 ? Float16(x.i) : Float16(Float32(x))
end
function Base.Float16(x::Normed{Ti,f}) where {Ti <: Union{UInt64, UInt128}, f}
    f == 1 ? Float16(x.i) : Float16(Float64(x))
end

function Base.Float32(x::Normed{UInt8,f}) where f
    f == 1 && return Float32(x.i)
    f == 2 && return Float32(Int32(x.i) * 0x101) * @f32(0x550055p-32)
    f == 3 && return Float32(Int32(x.i) * 0x00b) * @f32(0xd4c77bp-30)
    f == 4 && return Float32(Int32(x.i) * 0x101) * @f32(0x110011p-32)
    f == 5 && return Float32(Int32(x.i) * 0x003) * @f32(0xb02c0bp-30)
    f == 6 && return Float32(Int32(x.i) * 0x049) * @f32(0xe40039p-36)
    f == 7 && return Float32(Int32(x.i) * 0x01f) * @f32(0x852b5fp-35)
    f == 8 && return Float32(Int32(x.i) * 0x155) * @f32(0xc0f0fdp-40)
    0.0f0
end
function Base.Float32(x::Normed{UInt16,f}) where f
    f32 = Float32(x.i)
    f ==  1 && return f32
    f ==  2 && return f32 * @f32(0x55p-8)  + f32 * @f32(0x555555p-32)
    f ==  3 && return f32 * @f32(0x49p-9)  + f32 * @f32(0x249249p-33)
    f ==  4 && return f32 * @f32(0x11p-8)  + f32 * @f32(0x111111p-32)
    f ==  5 && return f32 * @f32(0x21p-10) + f32 * @f32(0x108421p-35)
    f ==  6 && return f32 * @f32(0x41p-12) + f32 * @f32(0x041041p-36)
    f ==  7 && return f32 * @f32(0x81p-14) + f32 * @f32(0x204081p-42)
    f == 16 && return f32 * @f32(0x01p-16) + f32 * @f32(0x010001p-48)
    Float32(x.i / rawone(x))
end
function Base.Float32(x::Normed{UInt32,f}) where f
    f == 1 && return Float32(x.i)
    i32 = unsafe_trunc(Int32, x.i)
    if f == 32
        rh, rl = Float32(i32>>>0x10), Float32((i32&0xFFFF)<<0x8 | i32>>>0x18)
        return muladd(rh, @f32(0x1p-16), rl * @f32(0x1p-40))
    elseif f >= 25
        rh, rl = Float32(i32>>>0x10), Float32((i32&0xFFFF)<<0xE + i32>>>UInt8(f-14))
        return muladd(rh, Float32(@exp2(16-f)), rl * Float32(@exp2(-14-f)))
    end
    # FIXME: avoid the branch in native x86_64 (non-SIMD) codes
    m = ifelse(i32 < 0, 0x1p32 * inv_rawone(x), 0.0)
    Float32(muladd(Float64(i32), inv_rawone(x), m))
end
function Base.Float32(x::Normed{Ti,f}) where {Ti <: Union{UInt64, UInt128}, f}
    f == 1 ? Float32(x.i) : Float32(Float64(x))
end

function Base.Float64(x::Normed{Ti,f}) where {Ti <: Union{UInt8, UInt16}, f}
    Float64(Normed{UInt32,f}(x))
end
function Base.Float64(x::Normed{UInt32,f}) where f
    f64 = Float64(x.i)
    f ==  1 && return f64
    f ==  2 && return (f64 * 0x040001p0) * 0x15555000015555p-72
    f ==  3 && return (f64 * 0x108421p0) * 0x11b6db76924929p-75
    f ==  4 && return (f64 * 0x010101p0) * 0x11000011000011p-72
    f ==  5 && return (f64 * 0x108421p0) * 0x04000002000001p-75
    f ==  6 && return (f64 * 0x09dfb1p0) * 0x1a56b8e38e6d91p-78
    f ==  7 && return (f64 * 0x000899p0) * 0x0f01480001e029p-70
    f ==  8 && return (f64 * 0x0a5a5bp0) * 0x18d300000018d3p-80
    f ==  9 && return (f64 * 0x001001p0) * 0x080381c8e3f201p-72
    f == 10 && return (f64 * 0x100001p0) * 0x04010000000401p-80
    f == 11 && return (f64 * 0x000009p0) * 0x0e3aaae3955639p-66
    f == 12 && return (f64 * 0x0a8055p0) * 0x186246e46e4cfdp-84
    f == 13 && return (f64 * 0x002001p0) * 0x10000004000001p-78
    f == 14 && return (f64 * 0x03400dp0) * 0x13b13b14ec4ec5p-84
    f == 15 && return (f64 * 0x000259p0) * 0x06d0c5a4f3a5e9p-75
    f == 16 && return (f64 * 0x011111p0) * 0x00f000ff00fff1p-80
    f == 18 && return (f64 * 0x0b06d1p0) * 0x17377445dd1231p-90
    f == 19 && return (f64 * 0x080001p0) * 0x00004000000001p-76
    f == 20 && return (f64 * 0x000101p0) * 0x0ff010ef10ff01p-80
    f == 21 && return (f64 * 0x004001p0) * 0x01fff8101fc001p-84
    f == 22 && return (f64 * 0x002945p0) * 0x18d0000000018dp-88
    f == 23 && return (f64 * 0x044819p0) * 0x07794a23729429p-92
    f == 27 && return (f64 * 0x000a21p0) * 0x0006518c7df9e1p-81
    f == 28 && return (f64 * 0x00000dp0) * 0x13b13b14ec4ec5p-84
    f == 30 && return (f64 * 0x001041p0) * 0x00fc003f03ffc1p-90
    f == 32 && return (f64 * 0x010101p0) * 0x00ff0000ffff01p-96
    f64 / rawone(x)
end
function Base.Float64(x::Normed{UInt64,f}) where f
    f == 1 && return Float64(x.i)
    if f >= 53
        rh = Float64(unsafe_trunc(Int64, x.i>>0x10)) * @exp2(16-f) # upper 48 bits
        rl = Float64(unsafe_trunc(Int32, x.i&0xFFFF)) * @exp2(-f)  # lower 16 bits
        return rh + muladd(rh, @exp2(-f), rl)
    end
    x.i / rawone(x)
end
function Base.Float64(x::Normed{UInt128,f}) where f
    f == 1 && return Float64(x.i)
    ih, il = unsafe_trunc(Int64, x.i>>0x40), unsafe_trunc(Int64, x.i)
    rh = Float64(ih>>>0x10) * @exp2(f <= 53 ? 80 : 80 - f) # upper 48 bits
    km = @exp2(f <= 53 ? 48 : 48 - f) # for middle 32 bits
    rm = Float64(unsafe_trunc(Int32, ih&0xFFFF)) * (0x1p16 * km) +
         Float64(unsafe_trunc(Int32, il>>>0x30)) * km
    rl = Float64(il&0xFFFFFFFFFFFF) * @exp2(f <= 53 ? 0 : -f) # lower 48 bits
    if f <= 53
        return (rh + (rm + rl)) / unsafe_trunc(Int64, rawone(x))
    elseif f < 76
        return rh + (rm + muladd(rh, @exp2(-f), rl))
    else
        return rh + (rm + rl)
    end
end

Base.BigFloat(x::Normed) = reinterpret(x) / BigFloat(rawone(x))

Base.Rational(x::Normed) = reinterpret(x)//rawone(x)

# unchecked arithmetic
*(x::T, y::T) where {T <: Normed} = convert(T,convert(floattype(T), x)*convert(floattype(T), y))
/(x::T, y::T) where {T <: Normed} = convert(T,convert(floattype(T), x)/convert(floattype(T), y))

# Functions
trunc(x::N) where {N <: Normed} = floor(x)
floor(x::N) where {N <: Normed} = reinterpret(N, x.i - x.i % rawone(N))
function ceil(x::Normed{T,f}) where {T, f}
    f == 1 && return x
    if typemax(T) % rawone(x) != 0
        upper = typemax(T) - typemax(T) % rawone(x)
        x.i > upper && throw_converterror(Normed{T,f}, ceil(T, typemax(x)))
    end
    r = x.i % rawone(x)
    reinterpret(Normed{T,f}, x.i - r + (r > 0 ? rawone(x) : zero(T)))
end
function round(x::Normed{T,f}) where {T, f}
    r = x.i % rawone(x)
    q = rawone(x) - r
    reinterpret(Normed{T,f}, r > q ? x.i + q : x.i - r)
end

trunc(::Type{Ti}, x::Normed) where {Ti <: Integer} = floor(Ti, x)
function floor(::Type{Ti}, x::Normed) where {Ti <: Integer}
    convert(Ti, reinterpret(x) ÷ rawone(x))
end
function ceil(::Type{Ti}, x::Normed) where {Ti <: Integer}
    d, r = divrem(x.i, rawone(x))
    convert(Ti, r > 0 ? d + oneunit(rawtype(x)) : d)
end
function round(::Type{Ti}, x::Normed) where {Ti <: Integer}
    d, r = divrem(x.i, rawone(x))
    convert(Ti, r > (rawone(x) >> 0x1) ? d + oneunit(rawtype(x)) : d)
end

# Iteration
# The main subtlety here is that iterating over N0f8(0):N0f8(1) will wrap around
# unless we iterate using a wider type
@inline start(r::StepRange{T}) where {T <: Normed} = widen1(reinterpret(r.start))
@inline next(r::StepRange{T}, i::Integer) where {T <: Normed} = (T(i,0), i+reinterpret(r.step))
@inline function done(r::StepRange{T}, i::Integer) where {T <: Normed}
    i1, i2 = reinterpret(r.start), reinterpret(r.stop)
    isempty(r) | (i < min(i1, i2)) | (i > max(i1, i2))
end

function decompose(x::Normed)
    g = gcd(reinterpret(x), rawone(x))
    div(reinterpret(x),g), 0, div(rawone(x),g)
end

# Range lengths
length(r::AbstractUnitRange{N}) where {N <: Normed{<:UShorterThanInt}} =
    floor(Int, last(r)) - floor(Int, first(r)) + 1
length(r::AbstractUnitRange{N}) where {N <: Normed{T}} where {T<:Unsigned} =
    r.start > r.stop ? T(0) : checked_add(floor(T, last(r)) - floor(T, first(r)), oneunit(T))

# Promotions
promote_rule(::Type{T}, ::Type{Tf}) where {T <: Normed,Tf <: AbstractFloat} = promote_type(floattype(T), Tf)
promote_rule(::Type{T}, ::Type{R}) where {T <: Normed,R <: Rational} = R
function promote_rule(::Type{T}, ::Type{Ti}) where {T <: Normed,Ti <: Union{Signed, Unsigned}}
    floattype(T)
end
@generated function promote_rule(::Type{Normed{T1,f1}}, ::Type{Normed{T2,f2}}) where {T1,T2,f1,f2}
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
    :(Normed{$T,$f})
end
