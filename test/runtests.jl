for f in ["ufixed.jl", "fixed32.jl"]
    println("Testing $f")
    include(f)
end
