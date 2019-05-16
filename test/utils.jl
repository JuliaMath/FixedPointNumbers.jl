""" return true if symbol `x` is a `ST` type"""
function _is_type(x::Symbol, ST::Symbol)
    try
        @eval $(x) isa Type && $(x) <: $(ST) && return true
    catch
        return false
    end
end

"""
    generate_fixedpoint_types(ST::Symbol=:FixedPoint)

generate a list of concrete `ST` types, where `ST ∈ (:FixedPoint, :Fixed, :Normed)`
"""
function generate_fixedpoint_types(ST::Symbol=:FixedPoint)
    fixed_types = setdiff(filter(x->_is_type(x, ST), names(FixedPointNumbers)), [:Fixed, :Normed, :FixedPoint])
    fixed_types = [@eval $(x) for x in fixed_types]
end
