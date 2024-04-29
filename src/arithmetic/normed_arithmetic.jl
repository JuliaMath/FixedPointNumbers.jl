using FixedPointNumbers: Normed
import .Wrapping: wrapping_mul
import .Saturating: saturating_mul
import .Checked: checked_mul

# wrapping arithmetic
function wrapping_mul(x::N, y::N) where {T <: Union{UInt8,UInt16,UInt32,UInt64}, f, N <: Normed{T,f}}
    z = widemul(x.i, y.i)
    N(div_2fm1(z, Val(Int(f))) % T, 0)
end

# saturating arithmetic
function saturating_mul(x::N, y::N) where {T <: Union{UInt8,UInt16,UInt32,UInt64}, f, N <: Normed{T,f}}
    f == bitwidth(T) && return wrapping_mul(x, y)
    z = min(widemul(x.i, y.i), widemul(typemax(N).i, rawone(N)))
    N(div_2fm1(z, Val(Int(f))) % T, 0)
end

# checked arithmetic
function checked_mul(x::N, y::N) where {N <: Normed}
    z = float(x) * float(y)
    z < typemax(N) + eps(N) / 2 || throw_overflowerror(:*, x, y)
    z % N
end
function checked_mul(x::N, y::N) where {T <: Union{UInt8,UInt16,UInt32,UInt64}, f, N <: Normed{T,f}}
    f == bitwidth(T) && return wrapping_mul(x, y)
    z = widemul(x.i, y.i)
    m = widemul(typemax(N).i, rawone(N)) + (rawone(N) >> 0x1)
    z < m || throw_overflowerror(:*, x, y)
    N(div_2fm1(z, Val(Int(f))) % T, 0)
end
