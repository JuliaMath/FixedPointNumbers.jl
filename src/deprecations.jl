import Base.@deprecate_binding

@deprecate_binding Fixed16 Q15f16
@deprecate_binding UFixed8 N0f8
@deprecate_binding UFixed10 N6f10
@deprecate_binding UFixed12 N4f12
@deprecate_binding UFixed14 N2f14
@deprecate_binding UFixed16 N0f16

@deprecate_binding UfixedBase Normed
@deprecate_binding Ufixed Normed
@deprecate_binding UFixed Normed
@deprecate_binding Ufixed8 N0f8
@deprecate_binding Ufixed10 N6f10
@deprecate_binding Ufixed12 N4f12
@deprecate_binding Ufixed14 N2f14
@deprecate_binding Ufixed16 N0f16

@deprecate_binding Fixed32 Q15f16

const UF = (N0f8, N6f10, N4f12, N2f14, N0f16)

@deprecate Fixed(x::Real) convert(Fixed{Int32, 16}, x)

@deprecate ufixed8(x)  N0f8(x)
@deprecate ufixed10(x) N6f10(x)
@deprecate ufixed12(x) N4f12(x)
@deprecate ufixed14(x) N2f14(x)
@deprecate ufixed16(x) N0f16(x)

## The next lines mimic the floating-point literal syntax "3.2f0"
# construction using a UInt, i.e., 0xccuf8
struct NormedConstructor{T,f} end
function *(n::Integer, ::NormedConstructor{T,f}) where {T,f}
    i = 8*sizeof(T)-f
    io = IOBuffer()
    show(io, n)
    nstr = String(take!(io))
    cstr = typeof(n) == T ? nstr : "convert($T, $nstr)"
    Base.depwarn("$(nstr)uf$f is deprecated, please use reinterpret(N$(i)f$f, $cstr) instead", :*)
    reinterpret(Normed{T,f}, convert(T, n))
end
const uf8  = NormedConstructor{UInt8,8}()
const uf10 = NormedConstructor{UInt16,10}()
const uf12 = NormedConstructor{UInt16,12}()
const uf14 = NormedConstructor{UInt16,14}()
const uf16 = NormedConstructor{UInt16,16}()

@deprecate_binding UfixedConstructor NormedConstructor
@deprecate_binding UFixedConstructor NormedConstructor
