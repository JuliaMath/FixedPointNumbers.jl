using FixedPointNumbers, Base.Test

@test isempty(detect_ambiguities(FixedPointNumbers, Base, Core))

for f in ["normed.jl", "fixed.jl"]
    println("Testing $f")
    include(f)
end
