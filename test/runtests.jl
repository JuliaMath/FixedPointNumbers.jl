using Aqua, FixedPointNumbers, Test

Aqua.test_all(FixedPointNumbers)

if Sys.ARCH === :x86_64 || Sys.ARCH === :i686
    using Documenter
    doctest(FixedPointNumbers, manual = false)
end

# older JET versions cannot analyze the package
if VERSION >= v"1.10"
    using JET
    @testset "JET" begin
        JET.test_package(FixedPointNumbers; target_modules=(FixedPointNumbers,))
    end
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
