using FixedPointNumbers, Test
using JET

@test isempty(detect_ambiguities(FixedPointNumbers, Base, Core))

@testset "JET" begin
    JET.test_package(FixedPointNumbers; target_modules=(FixedPointNumbers,))
end

@testset "normed" begin
    include("normed.jl")
end
@testset "fixed" begin
    include("fixed.jl")
end

@testset "traits" begin
    include("traits.jl")
end
