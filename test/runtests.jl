using FixedPointNumbers, Base.Test

for f in ["normed.jl", "fixed.jl"]
    println("Testing $f")
    include(f)
end

if v"0.5.0" <= VERSION < v"0.6.0-dev"
    @test isempty(detect_ambiguities(FixedPointNumbers, Base, Core))
end
