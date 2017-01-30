using BenchmarkTools
using FixedPointNumbers
import JLD

suite = BenchmarkGroup()
suite["conversion"] = BenchmarkGroup()

for FT in (Q0f7, Q1f14, Q7f24, N0f8, N2f14, N0f32, N2f30, N0f64, N2f62)
  x = FT(0.25)
  # Float16 doesn't behave well
  for T in (Float32, Float64, BigFloat)
    y = convert(T, x)
    suite["conversion"]["$FT-to-$T"] = @benchmarkable convert($T, $x)
    suite["conversion"]["$T-to-$FT"] = @benchmarkable convert($FT, $y)
  end
end

# Load the suite's cached parameters as part of including the file. This is much
# faster and more reliable than re-tuning `suite` every time the file is included
paramspath = Pkg.dir("FixedPointNumbers", "benchmark", "params.jld")
# tune!(suite); JLD.save(paramspath, "suite", params(suite));
loadparams!(suite, JLD.load(paramspath, "suite"), :evals, :samples);



