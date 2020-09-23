import Base.@deprecate_binding

function floattype(::Type{T}) where {T <: Real}
    Base.depwarn("""
        In a future release, the fallback definition of `floattype` will throw a MethodError if it cannot return a type `<:AbstractFloat`.
        See the documentation on `floattype` for guidance on whether to define a custom `floattype(::Type{$T})` method.
        """, :floattype)
    return T
end

function typechar(::Type{X}) where {X}
    Base.depwarn("""
        `typechar` was deprecated since the prefix may not be a single character in the future.
        We recommend not using private functions, but if you need to, use `type_prefix` instead.
        """, :typechar)
    Char(string(type_prefix(X))[1])
end
