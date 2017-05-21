using FixedPointNumbers, Base.Test

for f in ["normed.jl", "fixed.jl"]
    println("Testing $f")
    include(f)
end

@test isempty(detect_ambiguities(FixedPointNumbers, Base, Core))
