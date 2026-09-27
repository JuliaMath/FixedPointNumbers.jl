using FixedPointNumbers, Test
using FixedPointNumbers: bitwidth, rawtype

struct MyReal <: Real end

@testset "floattype" begin
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

    @test_skip(@test_throws MethodError floattype(MyReal))   # TODO: eliminate `@test_skipped` when depwarn is eliminated. See #177.
end

@testset "rawtype" begin
    @test rawtype(FixedPoint) === Union{}
    @test rawtype(Normed) === Union{}
    @test rawtype(Fixed) === Union{}
    @test rawtype(Normed{<:Unsigned}) === Union{}
    @test rawtype(Fixed{Union{Int8, Int16}}) === Union{Int8, Int16}
    @test rawtype(N0f8) === UInt8
    @test rawtype(Q0f7) === Int8
end
