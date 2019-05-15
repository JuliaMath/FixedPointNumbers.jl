@testset "floattype" begin
    function _is_fixed_type(x::Symbol)
        try
            @eval $(x) isa Type && $(x) <: FixedPoint && return true
        catch
            return false
        end
    end

    fixed_types = setdiff(filter(_is_fixed_type, names(FixedPointNumbers)), [:Fixed, :Normed, :FixedPoint])
    fixed_types = [@eval $(x) for x in fixed_types]

    exact_types = vcat([UInt8, UInt16, UInt32, UInt64, UInt128, Bool,
                       Int8, Int16, Int32, Int64, Int128],
                       fixed_types)
    for T in exact_types
        @test typemax(T) <= maxintfloat(floattype(T))
    end
end
