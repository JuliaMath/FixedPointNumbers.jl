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
