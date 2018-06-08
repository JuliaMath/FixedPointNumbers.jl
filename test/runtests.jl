using FixedPointNumbers, Test

@test isempty(detect_ambiguities(FixedPointNumbers, Base, Core))

@testset "normed" begin
    include("normed.jl")
end

@testset "fixed" begin
    include("fixed.jl")
end
