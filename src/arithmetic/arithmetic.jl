module FixedPointArithmetic

import Base: -, +, *, /, abs, div, fld, cld, rem, mod

import ..FixedPointNumbers
using FixedPointNumbers: FixedPoint
using FixedPointNumbers.Utilities

export mul_with_rounding

"""
    Wrapping

Submodule for wrapping arithmetic
"""
module Wrapping

import ..FixedPointNumbers
using FixedPointNumbers: FixedPoint
using FixedPointNumbers.Utilities

export wrapping_neg, wrapping_abs, wrapping_add, wrapping_sub, wrapping_mul
export wrapping_div, wrapping_fld, wrapping_cld, wrapping_rem, wrapping_mod
export wrapping_fdiv

wrapping_neg(x::X) where {X <: FixedPoint} = X(-x.i, 0)
wrapping_abs(x::X) where {X <: FixedPoint} = X(abs(x.i), 0)
wrapping_add(x::X, y::X) where {X <: FixedPoint} = X(x.i + y.i, 0)
wrapping_sub(x::X, y::X) where {X <: FixedPoint} = X(x.i - y.i, 0)
wrapping_mul(x::X, y::X) where {X <: FixedPoint} = (float(x) * float(y)) % X
function wrapping_fdiv(x::X, y::X) where {X <: FixedPoint}
    z = floattype(X)(x.i) / floattype(X)(y.i)
    isfinite(z) ? z % X : zero(X)
end
function wrapping_div(x::X, y::X, r::RoundingMode=RoundToZero) where {T,X <: FixedPoint{T}}
    z = round(floattype(X)(x.i) / floattype(X)(y.i), r)
    isfinite(z) || return zero(T)
    if T <: Unsigned
        _unsafe_trunc(T, z)
    else
        z > typemax(T) ? typemin(T) : _unsafe_trunc(T, z)
    end
end
wrapping_fld(x::X, y::X) where {X <: FixedPoint} = wrapping_div(x, y, RoundDown)
wrapping_cld(x::X, y::X) where {X <: FixedPoint} = wrapping_div(x, y, RoundUp)
wrapping_rem(x::X, y::X, r::RoundingMode=RoundToZero) where {T,X <: FixedPoint{T}} =
    X(x.i - wrapping_div(x, y, r) * y.i, 0)
wrapping_mod(x::X, y::X) where {X <: FixedPoint} = wrapping_rem(x, y, RoundDown)

end # module Wrapping

"""
    Saturating

Submodule for saturating arithmetic
"""
module Saturating

import ..FixedPointNumbers
using FixedPointNumbers: FixedPoint
using FixedPointNumbers.Utilities

export saturating_neg, saturating_abs, saturating_add, saturating_sub, saturating_mul
export saturating_div, saturating_fld, saturating_cld, saturating_rem, saturating_mod
export saturating_fdiv

saturating_neg(x::X) where {X <: FixedPoint} = X(~min(x.i - true, x.i), 0)
saturating_neg(x::X) where {X <: FixedPoint{<:Unsigned}} = zero(X)

saturating_abs(x::X) where {X <: FixedPoint} =
    X(ifelse(signbit(abs(x.i)), typemax(x.i), abs(x.i)), 0)

saturating_add(x::X, y::X) where {X <: FixedPoint} =
    X(x.i + ifelse(x.i < 0, max(y.i, typemin(x.i) - x.i), min(y.i, typemax(x.i) - x.i)), 0)
saturating_add(x::X, y::X) where {X <: FixedPoint{<:Unsigned}} = X(x.i + min(~x.i, y.i), 0)

saturating_sub(x::X, y::X) where {X <: FixedPoint} =
    X(x.i - ifelse(x.i < 0, min(y.i, x.i - typemin(x.i)), max(y.i, x.i - typemax(x.i))), 0)
saturating_sub(x::X, y::X) where {X <: FixedPoint{<:Unsigned}} = X(x.i - min(x.i, y.i), 0)

saturating_mul(x::X, y::X) where {X <: FixedPoint} = clamp(float(x) * float(y), X)

saturating_fdiv(x::X, y::X) where {X <: FixedPoint} =
    clamp(floattype(X)(x.i) / floattype(X)(y.i), X)

function saturating_div(x::X, y::X, r::RoundingMode=RoundToZero) where {T, X <: FixedPoint{T}}
    z = round(floattype(X)(x.i) / floattype(X)(y.i), r)
    isnan(z) && return zero(T)
    if T <: Unsigned
        isfinite(z) ? _unsafe_trunc(T, z) : typemax(T)
    else
        _unsafe_trunc(T, clamp(z, typemin(T), typemax(T)))
    end
end
saturating_fld(x::X, y::X) where {X <: FixedPoint} = saturating_div(x, y, RoundDown)
saturating_cld(x::X, y::X) where {X <: FixedPoint} = saturating_div(x, y, RoundUp)
function saturating_rem(x::X, y::X, r::RoundingMode=RoundToZero) where {T, X <: FixedPoint{T}}
    T <: Unsigned && r isa RoundingMode{:Up} && return zero(X)
    X(x.i - saturating_div(x, y, r) * y.i, 0)
end
saturating_mod(x::X, y::X) where {X <: FixedPoint} = saturating_rem(x, y, RoundDown)

end # module Saturating

"""
    CheckedArithMetic
"""
module Checked

import ..FixedPointNumbers
using FixedPointNumbers: FixedPoint, showtype
using FixedPointNumbers.Utilities

import Base.Checked: checked_neg, checked_abs, checked_add, checked_sub, checked_mul
import Base.Checked: checked_div, checked_fld, checked_cld, checked_rem, checked_mod

export checked_neg, checked_abs, checked_add, checked_sub, checked_mul
export checked_div, checked_fld, checked_cld, checked_rem, checked_mod
export checked_fdiv

checked_neg(x::X) where {X <: FixedPoint} = checked_sub(zero(X), x)
function checked_abs(x::X) where {X <: FixedPoint}
    abs(x.i) >= 0 || throw_overflowerror_abs(x)
    X(abs(x.i), 0)
end
function checked_add(x::X, y::X) where {X <: FixedPoint}
    r, f = Base.Checked.add_with_overflow(x.i, y.i)
    z = X(r, 0) # store first
    f && throw_overflowerror(:+, x, y)
    z
end
function checked_sub(x::X, y::X) where {X <: FixedPoint}
    r, f = Base.Checked.sub_with_overflow(x.i, y.i)
    z = X(r, 0) # store first
    f && throw_overflowerror(:-, x, y)
    z
end
function checked_mul(x::X, y::X) where {X <: FixedPoint}
    z = float(x) * float(y)
    typemin(X) - eps(X) / 2 <= z < typemax(X) + eps(X) / 2 || throw_overflowerror(:*, x, y)
    z % X
end
function checked_fdiv(x::X, y::X) where {T,X <: FixedPoint{T}}
    y === zero(X) && throw(DivideError())
    z = floattype(X)(x.i) / floattype(X)(y.i)
    if T <: Unsigned
        z < typemax(X) + eps(X) / 2 || throw_overflowerror(:/, x, y)
    else
        typemin(X) - eps(X) / 2 <= z < typemax(X) + eps(X) / 2 || throw_overflowerror(:/, x, y)
    end
    z % X
end
function checked_div(x::X, y::X, r::RoundingMode=RoundToZero) where {T, X <: FixedPoint{T}}
    y === zero(X) && throw(DivideError())
    z = round(floattype(X)(x.i) / floattype(X)(y.i), r)
    if T <: Signed
        z <= typemax(T) || throw_overflowerror_div(r, x, y)
    end
    _unsafe_trunc(T, z)
end
checked_fld(x::X, y::X) where {X <: FixedPoint} = checked_div(x, y, RoundDown)
checked_cld(x::X, y::X) where {X <: FixedPoint} = checked_div(x, y, RoundUp)
function checked_rem(x::X, y::X, r::RoundingMode=RoundToZero) where {T, X <: FixedPoint{T}}
    y === zero(X) && throw(DivideError())
    fx, fy = floattype(X)(x.i), floattype(X)(y.i)
    z = fx - round(fx / fy, r) * fy
    if T <: Unsigned && r isa RoundingMode{:Up}
        z >= zero(z) || throw_overflowerror_rem(r, x, y)
    end
    X(_unsafe_trunc(T, z), 0)
end
checked_mod(x::X, y::X) where {X <: FixedPoint} = checked_rem(x, y, RoundDown)

@noinline function throw_overflowerror(op::Symbol, @nospecialize(x), @nospecialize(y))
    io = IOBuffer()
    print(io, x, ' ', op, ' ', y, " overflowed for type ")
    showtype(io, typeof(x))
    throw(OverflowError(String(take!(io))))
end
@noinline function throw_overflowerror_abs(@nospecialize(x))
    io = IOBuffer()
    print(io, "abs(", x, ") overflowed for type ")
    showtype(io, typeof(x))
    throw(OverflowError(String(take!(io))))
end
@noinline function throw_overflowerror_div(r::RoundingMode, @nospecialize(x), @nospecialize(y))
    io = IOBuffer()
    op = r === RoundUp ? "cld(" : r === RoundDown ? "fld(" : "div("
    print(io, op, x, ", ", y, ") overflowed for type ", rawtype(x))
    throw(OverflowError(String(take!(io))))
end
@noinline function throw_overflowerror_rem(r::RoundingMode, @nospecialize(x), @nospecialize(y))
    io = IOBuffer()
    print(io, "rem(", x, ", ", y, ", ", r, ") overflowed for type ", typeof(x))
    throw(OverflowError(String(take!(io))))
end

end # module Checked

using .Wrapping
using .Saturating
using .Checked

using .Checked: throw_overflowerror

# re-export
for name in names(Wrapping)
    @eval export $name
end
for name in names(Saturating)
    @eval export $name
end
for name in names(Checked)
    @eval export $name
end

include("normed_arithmetic.jl")
include("fixed_arithmetic.jl")

# default arithmetic
const DEFAULT_ARITHMETIC = :wrapping

for (op, name) in ((:-, :neg), (:abs, :abs))
    f = Symbol(DEFAULT_ARITHMETIC, :_, name)
    @eval begin
        $op(x::X) where {X <: FixedPoint} = $f(x)
    end
end
for (op, name) in ((:+, :add), (:-, :sub), (:*, :mul))
    f = Symbol(DEFAULT_ARITHMETIC, :_, name)
    @eval begin
        $op(x::X, y::X) where {X <: FixedPoint} = $f(x, y)
    end
end

# force checked arithmetic
/(x::X, y::X) where {X <: FixedPoint} = checked_fdiv(x, y)
div(x::X, y::X, r::RoundingMode=RoundToZero) where {X <: FixedPoint} = checked_div(x, y, r)
fld(x::X, y::X) where {X <: FixedPoint} = checked_div(x, y, RoundDown)
cld(x::X, y::X) where {X <: FixedPoint} = checked_div(x, y, RoundUp)
rem(x::X, y::X) where {X <: FixedPoint} = checked_rem(x, y, RoundToZero)
rem(x::X, y::X, ::RoundingMode{:Down}) where {X <: FixedPoint} = checked_rem(x, y, RoundDown)
rem(x::X, y::X, ::RoundingMode{:Up}) where {X <: FixedPoint} = checked_rem(x, y, RoundUp)
mod(x::X, y::X) where {X <: FixedPoint} = checked_rem(x, y, RoundDown)

end # module FixedPointArithmetic
