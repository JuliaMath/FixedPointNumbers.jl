module FixedPointNumbersStatisticsExt

using FixedPointNumbers
import FixedPointNumbers: Treduce
using Statistics

import Statistics
if isdefined(Statistics, :_mean_promote)
    # https://github.com/JuliaMath/FixedPointNumbers.jl/pull/183
    Statistics._mean_promote(x::Real, y::FixedPoint) = Treduce(y)
end

end
