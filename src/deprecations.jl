import Base.@deprecate_binding

function floattype(::Type{T}) where {T <: Real}
    Base.depwarn("""
        In a future release, the fallback definition of `floattype` will throw a MethodError if it cannot return a type `<:AbstractFloat`.
        See the documentation on `floattype` for guidance on whether to define a custom `floattype(::Type{$T})` method.
        """, :floattype)
    return T
end
