using FixedPointNumbers: Fixed
import .Wrapping: wrapping_mul
import .Saturating: saturating_mul
import .Checked: checked_mul

# wrapping arithmetic
function wrapping_mul(x::F, y::F) where {T <: Union{Int8,Int16,Int32,Int64}, f, F <: Fixed{T,f}}
    z = widemul(x.i, y.i)
    F(div_2f(z, Val(Int(f))) % T, 0)
end

# saturating arithmetic
function saturating_mul(x::F, y::F) where {T <: Union{Int8,Int16,Int32,Int64}, f, F <: Fixed{T,f}}
    if f >= bitwidth(T) - 1
        z = min(widemul(x.i, y.i), widen1(typemax(T)) << f)
    else
        z = clamp(widemul(x.i, y.i), widen1(typemin(T)) << f, widen1(typemax(T)) << f)
    end
    F(div_2f(z, Val(Int(f))) % T, 0)
end

# checked arithmetic
function checked_mul(x::F, y::F) where {T <: Union{Int8,Int16,Int32,Int64}, f, F <: Fixed{T,f}}
    z = widemul(x.i, y.i)
    if f < 1
        m = widen1(typemax(T)) + 0x1
        n = widen1(typemin(T))
    else
        half = widen1(oneunit(T)) << (f - 1)
        m = widen1(typemax(T)) << f + half
        n = widen1(typemin(T)) << f - half
    end
    (n <= z) & (z < m) || throw_overflowerror(:*, x, y)
    F(div_2f(z, Val(Int(f))) % T, 0)
end

function mul_with_rounding(x::F, y::F, ::RoundingMode{:Nearest}) where {F <: Fixed}
    wrapping_mul(x, y)
end
function mul_with_rounding(x::F, y::F, ::RoundingMode{:NearestTiesUp}) where
{T<:Union{Int8,Int16,Int32,Int64},f,F <: Fixed{T,f}}
    z = widemul(x.i, y.i)
    F(((z + (oftype(z, 1) << f >>> 1)) >> f) % T, 0)
end
function mul_with_rounding(x::F, y::F, ::RoundingMode{:Down}) where
{T<:Union{Int8,Int16,Int32,Int64},f,F <: Fixed{T,f}}
    F((widemul(x.i, y.i) >> f) % T, 0)
end
