using FixedPointNumbers, Test

if VERSION >= v"1.6.0-DEV.816" # JuliaLang/julia #36962 # FIXME
    @test isempty(detect_ambiguities(FixedPointNumbers))
else
    @test isempty(detect_ambiguities(FixedPointNumbers, Base, Core))
end

using Documenter
doctest(FixedPointNumbers, manual = false)

@testset "normed" begin
    include("normed.jl")
end
@testset "fixed" begin
    include("fixed.jl")
end

@testset "traits" begin
    include("traits.jl")
end
