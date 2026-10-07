using FixedPointNumbers, Test
using FixedPointNumbers: bitwidth

struct MyReal <: Real end
struct MyInteger <: Integer end

@testset "floattype for non-concrete/non-FixedPoint types" begin
    @test floattype(FixedPoint) === BigFloat
    @test floattype(Normed) === BigFloat
    @test floattype(Fixed) === BigFloat
    @test floattype(FixedPoint{Int8}) === Float32
    @test floattype(Normed{UInt32}) === Float64
    @test floattype(Fixed{Int64}) === BigFloat
    @test floattype(FixedPoint{MyInteger}) === Float64
    @test floattype(FixedPoint{Union{},1}) === BigFloat
    @test floattype(Normed{Unsigned,2}) === BigFloat
    @test floattype(Fixed{T,3} where {T<:Signed}) === BigFloat
    @test_throws MethodError floattype(Fixed{T,3} where T) # T is not bounded

    for T in (UInt8, UInt16, UInt32, UInt64, UInt128, Bool,
              Int8, Int16, Int32, Int64, Int128)
        @test typemax(T) <= maxintfloat(floattype(T))
    end
    @test floattype(Rational{Int}) === Float64
    @test floattype(Complex{Int16})   === Complex{Float32}
    @test floattype(Complex{Float32}) === Complex{Float32}
    @test floattype(Base.TwicePrecision{Float16}) === Float32
    @test floattype(Base.TwicePrecision{Float32}) === Float64
    @test floattype(Base.TwicePrecision{Float64}) === Float64
    @test floattype(typeof(π))                    === Float64

    @test_throws MethodError floattype(MyReal) # See #177.
    @test floattype(MyInteger) === Float64
end
