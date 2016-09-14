import Base.@deprecate_binding

@deprecate_binding UfixedBase UFixed
@deprecate_binding UfixedConstructor UFixedConstructor
@deprecate_binding Ufixed UFixed
@deprecate_binding Ufixed8 UFixed8
@deprecate_binding Ufixed10 UFixed10
@deprecate_binding Ufixed12 UFixed12
@deprecate_binding Ufixed14 UFixed14
@deprecate_binding Ufixed16 UFixed16

@deprecate_binding Fixed32 Fixed16
@deprecate Fixed(x::Real) convert(Fixed{Int32, 16}, x)

@deprecate ufixed8(x)  UFixed8(x)
@deprecate ufixed10(x) UFixed10(x)
@deprecate ufixed12(x) UFixed12(x)
@deprecate ufixed14(x) UFixed14(x)
@deprecate ufixed16(x) UFixed16(x)

Compat.@dep_vectorize_1arg Real ufixed8
Compat.@dep_vectorize_1arg Real ufixed10
Compat.@dep_vectorize_1arg Real ufixed12
Compat.@dep_vectorize_1arg Real ufixed14
Compat.@dep_vectorize_1arg Real ufixed16
