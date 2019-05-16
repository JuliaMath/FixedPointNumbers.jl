@testset "floattype" begin
    exact_types = vcat([UInt8, UInt16, UInt32, UInt64, UInt128, Bool,
                       Int8, Int16, Int32, Int64, Int128],
                       generate_fixedpoint_types())
    for T in exact_types
        @test typemax(T) <= maxintfloat(floattype(T))
    end
end
