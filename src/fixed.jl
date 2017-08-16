# 32-bit fixed point; parameter `f` is the number of fraction bits
struct Fixed{T <: Signed,f} <: FixedPoint{T,  f}
    i::T

    # constructor for manipulating the representation;
    # selected by passing an extra dummy argument
    (::Type{Fixed{T, f}})(i::Integer, _) where {T,f} = new{T, f}(i % T)
    (::Type{Fixed{T, f}})(x) where {T,f} = convert(Fixed{T,f}, x)
end

reinterpret(::Type{Fixed{T,f}}, x::T) where {T <: Signed,f} = Fixed{T,f}(x, 0)

typechar(::Type{X}) where {X <: Fixed} = 'Q'
signbits(::Type{X}) where {X <: Fixed} = 1

for T in (Int8, Int16, Int32, Int64)
    for f in 0:sizeof(T)*8-1
        sym = Symbol(String(take!(showtype(_iotypealias, Fixed{T,f}))))
        @eval begin
            const $sym = Fixed{$T,$f}
            export $sym
        end
    end
end

# basic operators
-(x::Fixed{T,f}) where {T,f} = Fixed{T,f}(-x.i,0)
abs(x::Fixed{T,f}) where {T,f} = Fixed{T,f}(abs(x.i),0)

+(x::Fixed{T,f}, y::Fixed{T,f}) where {T,f} = Fixed{T,f}(x.i+y.i,0)
-(x::Fixed{T,f}, y::Fixed{T,f}) where {T,f} = Fixed{T,f}(x.i-y.i,0)

# with truncation:
#*{f}(x::Fixed32{f}, y::Fixed32{f}) = Fixed32{f}(Base.widemul(x.i,y.i)>>f,0)
# with rounding up:
*(x::Fixed{T,f}, y::Fixed{T,f}) where {T,f} = Fixed{T,f}((Base.widemul(x.i,y.i) + (convert(widen(T), 1) << (f-1) ))>>f,0)

/(x::Fixed{T,f}, y::Fixed{T,f}) where {T,f} = Fixed{T,f}(div(convert(widen(T), x.i) << f, y.i), 0)


# # conversions and promotions
convert(::Type{Fixed{T,f}}, x::Integer) where {T,f} = Fixed{T,f}(round(T, convert(widen1(T),x)<<f),0)
convert(::Type{Fixed{T,f}}, x::AbstractFloat) where {T,f} = Fixed{T,f}(round(T, trunc(widen1(T),x)<<f + rem(x,1)*(1<<f)),0)
convert(::Type{Fixed{T,f}}, x::Rational) where {T,f} = Fixed{T,f}(x.num)/Fixed{T,f}(x.den)

rem(x::Integer, ::Type{Fixed{T,f}}) where {T,f} = Fixed{T,f}(rem(x,T)<<f,0)
rem(x::Real,    ::Type{Fixed{T,f}}) where {T,f} = Fixed{T,f}(rem(Integer(trunc(x)),T)<<f + rem(Integer(round(rem(x,1)*(1<<f))),T),0)

# convert{T,f}(::Type{AbstractFloat}, x::Fixed{T,f}) = convert(floattype(x), x)
float(x::Fixed) = convert(floattype(x), x)

convert(::Type{BigFloat}, x::Fixed{T,f}) where {T,f} =
    convert(BigFloat,x.i>>f) + convert(BigFloat,x.i&(1<<f - 1))/convert(BigFloat,1<<f)
convert(::Type{TF}, x::Fixed{T,f}) where {TF <: AbstractFloat,T,f} =
    convert(TF,x.i>>f) + convert(TF,x.i&(1<<f - 1))/convert(TF,1<<f)

convert(::Type{Bool}, x::Fixed{T,f}) where {T,f} = x.i!=0
function convert(::Type{Integer}, x::Fixed{T,f}) where {T,f}
    isinteger(x) || throw(InexactError())
    convert(Integer, x.i>>f)
end
function convert(::Type{TI}, x::Fixed{T,f}) where {TI <: Integer,T,f}
    isinteger(x) || throw(InexactError())
    convert(TI, x.i>>f)
end

convert(::Type{TR}, x::Fixed{T,f}) where {TR <: Rational,T,f} =
    convert(TR, x.i>>f + (x.i&(1<<f-1))//(1<<f))

promote_rule(ft::Type{Fixed{T,f}}, ::Type{TI}) where {T,f,TI <: Integer} = Fixed{T,f}
promote_rule(::Type{Fixed{T,f}}, ::Type{TF}) where {T,f,TF <: AbstractFloat} = TF
promote_rule(::Type{Fixed{T,f}}, ::Type{Rational{TR}}) where {T,f,TR} = Rational{TR}

# TODO: Document and check that it still does the right thing.
decompose(x::Fixed{T,f}) where {T,f} = x.i, -f, 1
