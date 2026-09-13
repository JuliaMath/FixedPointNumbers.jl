module FixedPointArithmetic

using ..FixedPointNumbers

import Base: -, +, *, /, abs, div, fld, cld, rem, mod

module Wrapping

export wrapping_neg, wrapping_abs, wrapping_add, wrapping_sub, wrapping_mul,
       wrapping_div, wrapping_fld, wrapping_cld, wrapping_rem, wrapping_mod,
       wrapping_fdiv

for name in names(Wrapping)
    startswith(string(name), "wrapping_") || continue
    @eval function $name end
end

end # module Wrapping

module Saturating

export saturating_neg, saturating_abs, saturating_add, saturating_sub, saturating_mul,
       saturating_div, saturating_fld, saturating_cld, saturating_rem, saturating_mod,
       saturating_fdiv

for name in names(Saturating)
    startswith(string(name), "saturating_") || continue
    @eval function $name end
end

end # module Saturating

module Checked

import Base.Checked: checked_neg, checked_abs, checked_add, checked_sub, checked_mul,
                     checked_div, checked_fld, checked_cld, checked_rem, checked_mod

export checked_neg, checked_abs, checked_add, checked_sub, checked_mul,
       checked_div, checked_fld, checked_cld, checked_rem, checked_mod,
       checked_fdiv

function checked_fdiv end

end # module Checked

module Unchecked

using ..FixedPointNumbers
using ..Wrapping

export unchecked_neg, unchecked_abs, unchecked_add, unchecked_sub, unchecked_mul,
       unchecked_div, unchecked_fld, unchecked_cld, unchecked_rem, unchecked_mod,
       unchecked_fdiv

for name in (:neg, :abs)
    fu = Symbol(:unchecked_, name)
    fw = Symbol(:wrapping_, name)
    @eval begin
        $fu(x::X) where {X <: FixedPoint} = $fw(x)
    end
end
for name in (:add, :sub, :mul, :div, :fld, :cld, :rem, :mod, :fdiv)
    fu = Symbol(:unchecked_, name)
    fw = Symbol(:wrapping_, name)
    @eval begin
        $fu(x::X, y::X) where {X <: FixedPoint} = $fw(x, y)
    end
    name in (:div, :rem) || continue
    @eval begin
        $fu(x::X, y::X, r::RoundingMode{M}) where {X <: FixedPoint, M} = $fw(x, y, r)
    end
end

end # module Unchecked

using .Wrapping, .Saturating, .Checked, .Unchecked

# re-export
for Mod in (Wrapping, Saturating, Checked, Unchecked)
    for name in names(Mod)
        @eval export $name
    end
end

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

const DEFAULT_DIV_ARITHMETIC = :checked

for name in (:fdiv, :div, :fld, :cld, :rem, :mod)
    f = Symbol(DEFAULT_DIV_ARITHMETIC, :_, name)
    if name === :fdiv
        @eval begin
            /(x::X, y::X) where {X <: FixedPoint} = $f(x, y)
        end
        continue
    end
    @eval begin
        $name(x::X, y::X) where {X <: FixedPoint} = $f(x, y)
    end
    name in (:div, :rem) || continue
end

for m in (:(:Nearest), :(:ToZero), :(:Up), :(:Down))
    _div = Symbol(DEFAULT_DIV_ARITHMETIC, :_div)
    _rem = Symbol(DEFAULT_DIV_ARITHMETIC, :_rem)
    @eval begin
        div(x::X, y::X, r::RoundingMode{$m}) where {X <: FixedPoint} = $_div(x, y, r)
        rem(x::X, y::X, r::RoundingMode{$m}) where {X <: FixedPoint} = $_rem(x, y, r)
    end
end

end # module FixedPointArithmetic
