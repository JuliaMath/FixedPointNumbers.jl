using FixedPointNumbers, Base.Test

if VERSION >= v"0.5.0"
    @test isempty(detect_ambiguities(FixedPointNumbers, Base, Core))
end

for f in ["normed.jl", "fixed.jl"]
    println("Testing $f")
    include(f)
end
