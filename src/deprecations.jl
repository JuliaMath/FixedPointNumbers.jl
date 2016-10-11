import Base.@deprecate_binding

@deprecate_binding Fixed16 Q15f16
@deprecate_binding UFixed8 N0f8
@deprecate_binding UFixed10 N6f10
@deprecate_binding UFixed12 N4f12
@deprecate_binding UFixed14 N2f14
@deprecate_binding UFixed16 N0f16

@deprecate_binding UfixedBase UFixed
@deprecate_binding UfixedConstructor UFixedConstructor
@deprecate_binding Ufixed UFixed
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

Compat.@dep_vectorize_1arg Real ufixed8
Compat.@dep_vectorize_1arg Real ufixed10
Compat.@dep_vectorize_1arg Real ufixed12
Compat.@dep_vectorize_1arg Real ufixed14
Compat.@dep_vectorize_1arg Real ufixed16
