using FixedPointNumbers, Test

@test isempty(detect_ambiguities(FixedPointNumbers, Base, Core))

if Sys.ARCH === :x86_64 || Sys.ARCH === :i686
    using Documenter
    doctest(FixedPointNumbers, manual = false)
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
