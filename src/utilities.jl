# utility functions and macros, which are independent of `FixedPoint`
bitwidth(T::Type) = 8sizeof(T)

widen1(::Type{Int8})    = Int16
widen1(::Type{UInt8})   = UInt16
widen1(::Type{Int16})   = Int32
widen1(::Type{UInt16})  = UInt32
widen1(::Type{Int32})   = Int64
widen1(::Type{UInt32})  = UInt64
widen1(::Type{Int64})   = Int128
widen1(::Type{UInt64})  = UInt128
widen1(::Type{Int128})  = Int128
widen1(::Type{UInt128}) = UInt128
widen1(x::Integer) = x % widen1(typeof(x))

const ShortInts = Union{Int8, UInt8, Int16, UInt16}
const LongInts = Union{Int64, UInt64, Int128, UInt128, BigInt}

const ShorterThanInt = Int === Int32 ? ShortInts : Union{ShortInts, Int32, UInt32}
const NotBiggerThanInt = Union{ShorterThanInt, Int, UInt}
const SShorterThanInt = typeintersect(ShorterThanInt, Signed)
const UShorterThanInt = typeintersect(ShorterThanInt, Unsigned)

macro f32(x::Float64) # just for hexadecimal floating-point literals
    :(Float32($x))
end
macro exp2(n)
     :(_exp2(Val($(esc(n)))))
end
_exp2(::Val{N}) where {N} = exp2(N)

# these are defined in julia/float.jl or julia/math.jl, but not exported
significand_bits(::Type{Float32}) = 23
significand_bits(::Type{Float64}) = 52
exponent_bias(::Type{Float32}) = 127
exponent_bias(::Type{Float64}) = 1023
