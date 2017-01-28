using FixedPointNumbers, Base.Test

for f in ["normed.jl", "fixed.jl"]
    println("Testing $f")
    include(f)
end

if VERSION >= v"0.5.0"
    @test isempty(detect_ambiguities(FixedPointNumbers, Base, Core))
end
