import Base.@deprecate_binding

function typechar(::Type{X}) where {X}
    Base.depwarn("""
        `typechar` was deprecated since the prefix may not be a single character in the future.
        We recommend not using private functions, but if you need to, use `type_prefix` instead.
        """, :typechar)
    Char(string(type_prefix(X))[1])
end
