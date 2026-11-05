using FixedPointNumbers, Test
using FixedPointNumbers: bitwidth

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

@testset "`Union{}` as type argument" begin
    @test_throws MethodError floattype(Union{})
    @test_throws MethodError scaledual(Union{}, 0.5N0f8)
    @test_throws MethodError 0.5 * Union{}
    @test_throws MethodError 0.5 % Union{}
    for f in (wrapping_mul, saturating_mul, Base.Checked.checked_mul,
              wrapping_rem, saturating_rem, Base.Checked.checked_rem)
        @test_throws MethodError f(0.5, Union{})
    end
end
