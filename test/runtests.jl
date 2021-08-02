using FixedPointNumbers, Test

if VERSION >= v"1.6.0-DEV.816" # JuliaLang/julia #36962 # FIXME
    @test isempty(detect_ambiguities(FixedPointNumbers))
else
    @test isempty(detect_ambiguities(FixedPointNumbers, Base, Core))
end

if Sys.ARCH === :x86_64 || Sys.ARCH === :i686
    using Documenter
    doctest(FixedPointNumbers, manual = false)
end

if !isdefined(Base, :checked_length)
    const checked_length = Base.unsafe_length
else
    using Base: checked_length
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
