using FixedPointNumbers, Statistics, Test
using FixedPointNumbers: bitwidth

@testset "domain of f" begin
    @test_throws DomainError zero(Normed{UInt8,-1})
    @test_throws DomainError zero(Normed{UInt8,0})
    @test_throws DomainError zero(Normed{UInt8,9})
    @test_throws DomainError zero(Normed{UInt16,17})
end

@testset "reinterpret/bitstring" begin
    @test reinterpret(N0f8, 0xa2).i  === 0xa2
    @test reinterpret(N6f10, 0x1fa2).i === 0x1fa2
    @test reinterpret(N4f12, 0x1fa2).i === 0x1fa2
    @test reinterpret(N2f14, 0x1fa2).i === 0x1fa2
    @test reinterpret(N0f16, 0x1fa2).i === 0x1fa2

    @test reinterpret(reinterpret(N0f8, 0xa2))    === 0xa2
    @test reinterpret(reinterpret(N6f10, 0x00a2)) === 0x00a2
    @test reinterpret(reinterpret(N4f12, 0x00a2)) === 0x00a2
    @test reinterpret(reinterpret(N2f14, 0x00a2)) === 0x00a2
    @test reinterpret(reinterpret(N0f16, 0x00a2)) === 0x00a2

    @test reinterpret(UInt8, 1N0f8) === 0xff

    @test bitstring(reinterpret(N0f8, 0xa2))    === "10100010"
    @test bitstring(reinterpret(N6f10, 0x00a2)) === "0000000010100010"

    @test 0.635N0f8   == N0f8(0.635)
    @test 0.635N6f10 == N6f10(0.635)
    @test 0.635N4f12 == N4f12(0.635)
    @test 0.635N2f14 == N2f14(0.635)
    @test 0.635N0f16 == N0f16(0.635)

    @test N0f8(1.0) == reinterpret(N0f8, 0xff)
    @test N0f8(0.5) == reinterpret(N0f8, 0x80)
    @test N2f14(1.0) == reinterpret(N2f14, 0x3fff)
    v = N4f12.([2])
    @test v == N4f12[reinterpret(N4f12, 0x1ffe)]
    @test isa(v, Vector{N4f12})
end

UF2 = (Normed{UInt32,16}, Normed{UInt64,3}, Normed{UInt64,51}, Normed{UInt128,7}, Normed{UInt128,51})

@testset "limits and identities" begin
    for T in (FixedPointNumbers.UF..., UF2...)
        @test zero(T) == 0
        @test one(T) == 1
        @test one(T) * one(T) == one(T)
        @test typemin(T) == 0
        @test floatmin(T) == eps(T)
        @test eps(zero(T)) == eps(typemax(T))
        @test sizeof(T) == sizeof(FixedPointNumbers.rawtype(T))
    end
    @test typemax(N0f8) == 1
    @test typemax(N6f10) == typemax(UInt16)//(2^10-1)
    @test typemax(N4f12) == typemax(UInt16)//(2^12-1)
    @test typemax(N2f14) == typemax(UInt16)//(2^14-1)
    @test typemax(N0f16) == 1
    @test typemax(N6f10) == typemax(UInt16) // (2^10-1)
    @test typemax(N4f12) == typemax(UInt16) // (2^12-1)
    @test typemax(N2f14) == typemax(UInt16) // (2^14-1)
    @test typemax(Normed{UInt32,16}) == typemax(UInt32) // (2^16-1)
    @test typemax(Normed{UInt64,3}) == typemax(UInt64) // (2^3-1)
    @test typemax(Normed{UInt128,7}) == typemax(UInt128) // (2^7-1)
    @test typemax(Normed{UInt128,100}) == typemax(UInt128) // (UInt128(2)^100-1)
end

@testset "inexactness" begin
    # TODO: change back to InexactError when it allows message strings
    @test_throws ArgumentError N0f8(2)
    @test_throws ArgumentError N0f8(255)
    @test_throws ArgumentError N0f8(0xff)
    @test_throws ArgumentError N0f16(2)
    @test_throws ArgumentError N0f16(0xff)
    @test_throws ArgumentError N0f16(0xffff)
    @test_throws ArgumentError convert(N0f8,  typemax(N6f10))
    @test_throws ArgumentError convert(N0f16, typemax(N6f10))
    @test_throws ArgumentError convert(Normed{UInt128,100}, 10^9)
    @test_throws ArgumentError convert(Normed{UInt128,100}, 10.0^9)
end

@testset "conversion" begin
    x = N0f8(0.5)
    @test convert(N0f8, x) === x

    @test convert(N0f8,  1.1/typemax(UInt8)) == eps(N0f8)
    @test convert(N6f10, 1.1/typemax(UInt16)*64) == eps(N6f10)
    @test convert(N4f12, 1.1/typemax(UInt16)*16) == eps(N4f12)
    @test convert(N2f14, 1.1/typemax(UInt16)*4)  == eps(N2f14)
    @test convert(N0f16, 1.1/typemax(UInt16))    == eps(N0f16)
    @test convert(Normed{UInt32,16}, 1.1/typemax(UInt32)*2^16) == eps(Normed{UInt32,16})
    @test convert(Normed{UInt64,3},  1.1/typemax(UInt64)*UInt64(2)^61)  == eps(Normed{UInt64,3})
    @test convert(Normed{UInt128,7}, 1.1/typemax(UInt128)*UInt128(2)^121) == eps(Normed{UInt128,7})

    @test convert(N0f8,  1.1f0/typemax(UInt8)) == eps(N0f8)

    @test convert(N0f8, 1//255) === eps(N0f8)
    @test convert(N0f8, Rational{Int8}(3//5)) === N0f8(3/5)
    @test convert(N0f8, Rational{UInt8}(3//5)) === N0f8(3/5)
    @test_throws ArgumentError convert(N0f8, typemax(Rational{UInt8}))

    @test convert(N0f8, Base.TwicePrecision(1.0)) === 1N0f8

    @test convert(Float64, eps(N0f8)) == 1/typemax(UInt8)
    @test convert(Float32, eps(N0f8)) == 1.0f0/typemax(UInt8)
    @test convert(BigFloat, eps(N0f8)) == BigFloat(1)/typemax(UInt8)
    for T in (FixedPointNumbers.UF..., UF2...)
        @test convert(Bool, zero(T)) == false
        @test convert(Bool, one(T))  == true
        @test_throws InexactError convert(Bool, convert(T, 0.2))
        @test convert(Int, one(T)) == 1
        @test convert(Integer, one(T)) == 1
        @test convert(Rational, one(T)) == 1
    end
    @test convert(N0f16, one(N0f8)) === one(N0f16)
    @test convert(N0f16, N0f8(0.5)).i === 0x8080
    @test convert(Normed{UInt16,7}, Normed{UInt8,7}(0.504)) === Normed{UInt16,7}(0.504)
end

@testset "integer conversions" begin
    @test convert(UInt, 1N1f7) === UInt(1)
    @test convert(Integer, 1N1f7) === 0x01
    @test convert(Int, 1N1f7) === 1
    @test_throws InexactError convert(Integer, 0.5N1f7)
    @test_throws InexactError convert(Int8, 256N8f8)
end

@testset "rational conversions" begin
    @test convert(Rational, 0.5N0f8) === Rational{UInt8}(0x80//0xff)
    @test convert(Rational, 0.5N4f12) === Rational{UInt16}(0x800//0xfff)
    @test convert(Rational{Int}, 0.5N0f8) === Rational{Int}(0x80//0xff)

    @test rationalize(0.8N0f8) === Rational{Int}(4//5)
    @test rationalize(Int16, 0.804N0f8) === Rational{Int16}(41//51)
    @test rationalize(0.804N0f8, tol=0.002) === Rational{Int}(41//51)
    @test rationalize(Int8, 0.804N0f8, tol=0.005) === Rational{Int8}(4//5)
end

@testset "BigFloat conversions" begin
    @test convert(BigFloat, 0.5N0f8)::BigFloat == 128 / big"255"

    @test big(N7f1) === BigFloat # !== BigInt
    @test big(0.5N4f4)::BigFloat == 8 / big"15"
end

@testset "conversion from float" begin
    # issue 102
    for T in (UInt8, UInt16, UInt32, UInt64, UInt128)
        for Tf in (Float16, Float32, Float64)
            @testset "Normed{$T,$f}(::$Tf)" for f = 1:bitwidth(T)
                N = Normed{T,f}
                r = FixedPointNumbers.rawone(N)

                @test reinterpret(N(zero(Tf))) == 0x0

                input_typemax = Tf(typemax(N))
                if isinf(input_typemax)
                    @test reinterpret(N(floatmax(Tf))) >= round(T, floatmax(Tf))
                else
                    @test reinterpret(N(input_typemax)) > (typemax(T)>>1) # overflow check
                    @test N(input_typemax) >= N(prevfloat(input_typemax))
                end

                input_upper = Tf(BigFloat(typemax(T)) / r, RoundDown)
                isinf(input_upper) && continue # for Julia v0.7
                @test reinterpret(N(input_upper)) == T(min(round(BigFloat(input_upper) * r), typemax(T)))

                input_exp2 = Tf(exp2(bitwidth(T) - f))
                isinf(input_exp2) && continue
                @test reinterpret(N(input_exp2)) == T(input_exp2) * r
            end
        end
    end
    @test N0f32(Float32(0x0.7FFFFFp-32)) == zero(N0f32)
    @test N0f32(Float32(0x0.800000p-32)) <= eps(N0f32) # should be zero in RoundNearest mode
    @test N0f32(Float32(0x0.800001p-32)) == eps(N0f32)
end

@testset "conversions to float" begin
    x = N0f8(0.3)
    for T in (Float16, Float32, Float64)
        y = convert(T, x)
        @test isa(y, T)
    end

    for Tf in (Float16, Float32, Float64)
        @testset "$Tf(::Normed{$T})" for T in (UInt8, UInt16)
            @testset "$Tf(::Normed{$T,$f})" for f = 1:bitwidth(T)
                N = Normed{T,f}
                float_err = 0.0
                for i = typemin(T):typemax(T)
                    f_expected = Tf(i / BigFloat(FixedPointNumbers.rawone(N)))
                    isinf(f_expected) && break # for Float16(::Normed{UInt16,1})
                    f_actual = Tf(reinterpret(N, i))
                    float_err += abs(f_actual - f_expected)
                end
                @test float_err == 0.0
            end
        end
        @testset "$Tf(::Normed{$T})" for T in (UInt32, UInt64, UInt128)
            @testset "$Tf(::Normed{$T,$f})" for f = 1:bitwidth(T)
                N = Normed{T,f}
                error_count = 0
                for i in vcat(T(0x00):T(0xFF), (typemax(T)-0xFF):typemax(T))
                    f_expected = Tf(i / BigFloat(FixedPointNumbers.rawone(N)))
                    isinf(f_expected) && break # for Float16() and Float32()
                    f_actual = Tf(reinterpret(N, i))
                    f_actual == f_expected && continue
                    f_actual == prevfloat(f_expected) && continue
                    f_actual == nextfloat(f_expected) && continue
                    error_count += 1
                end
                @test error_count == 0
            end
        end
    end
end

@testset "modulus" begin
    @test  N0f8(0.2) % N0f8  === N0f8(0.2)
    @test N2f14(1.2) % N0f16 === N0f16(0.20002)
    @test N2f14(1.2) % N0f8  === N0f8(0.196)

    for i = 0.0:0.1:1.0
        @test i % N0f8 === N0f8(i)
    end
    @test ( 1.5 % N0f8).i == round(Int,  1.5*255) % UInt8
    @test (-0.3 % N0f8).i == round(Int, -0.3*255) % UInt8

    for i = 0.0:0.1:64.0
        @test i % N6f10 === N6f10(i)
    end
    @test (65.2 % N6f10).i == round(Int, 65.2*1023) % UInt16
    @test (-0.3 % N6f10).i == round(Int, -0.3*1023) % UInt16

    @test 1 % N0f8 == 1
    @test 2 % N0f8 == N0f8(0.996)

    # issue #150
    @test all(f -> 1.0f0 % Normed{UInt32,f} == oneunit(Normed{UInt32,f}), 1:32)
    @test all(f -> 1.0e0 % Normed{UInt64,f} == oneunit(Normed{UInt64,f}), 1:64)
end

@testset "bitwise" begin
    x = N0f8(0b01010001, 0)
    @test ~x == N0f8(0b10101110, 0)
    @test -x == reinterpret(N0f8, 0xaf)
end

@testset "float" begin
    @test isa(float(one(Normed{UInt8,7})),   Float32)
    @test isa(float(one(Normed{UInt32,18})), Float64)
    @test isa(float(one(Normed{UInt32,25})), Float64)
end

@testset "arithmetic" begin
    for T in (FixedPointNumbers.UF..., UF2...)
        x = T(0x10,0)
        y = T(0x25,0)
        fx = float(x)
        fy = float(y)
        @test y > x
        @test y != x
        @test typeof(x+y) == T
        @test typeof((x+y)-y) == T
        @test typeof(x*y) == T
        @test typeof(x/y) == T
        @test (x+y) ≈ T(0x35,0)
        @test ((x+y)-x) ≈ fy
        @test ((x-y)+y) ≈ fx
        @test (x*y)  ≈ convert(T, fx*fy)
        @test (x/y)  ≈ convert(T, fx/fy)
        @test (x^2) ≈ convert(T, fx^2)
        @test (x^2.1f0) ≈ fx^2.1f0
        @test (x^2.1) ≈ convert(Float64, x)^2.1
    end
end

@testset "rounding" begin
    for T in (UInt8, UInt16, UInt32, UInt64)
        rs = vcat([ oneunit(T) << b - oneunit(T) << 1 for b = 1:bitwidth(T)],
                  [ oneunit(T) << b - oneunit(T)      for b = 1:bitwidth(T)],
                  [ oneunit(T) << b                   for b = 2:bitwidth(T)-1])
        @testset "rounding Normed{$T,$f}" for f = 1:bitwidth(T)
            N = Normed{T,f}
            xs = (reinterpret(N, r) for r in rs)
            @test all(x -> trunc(x) == trunc(float(x)), xs)
            @test all(x -> floor(x) == floor(float(x)), xs)
            # force `Normed` comparison avoiding rounding errors
            @test all(x -> ceil(float(x)) > typemax(N) || ceil(x) == N(ceil(float(x))), xs)
            @test all(x -> round(x) == round(float(x)), xs)
            @test all(x -> trunc(UInt64, x) === trunc(UInt64, float(x)), xs)
            @test all(x -> floor(UInt64, x) === floor(UInt64, float(x)), xs)
            @test all(x ->  ceil(UInt64, x) ===  ceil(UInt64, float(x)), xs)
            @test all(x -> round(UInt64, x) === round(UInt64, float(x)), xs)
        end
    end
    @testset "rounding Normed{Int16,$f} with overflow" for f in filter(x->!ispow2(x), 1:16)
        N = Normed{UInt16,f}
        @test_throws ArgumentError ceil(typemax(N))
        @test_throws ArgumentError ceil(floor(typemax(N)) + eps(N))
    end
    @testset "rounding mode" begin
        @test round(1.504N1f7, RoundNearest) === 2N1f7
        @test round(1.504N1f7, RoundToZero) === 1N1f7
        @test round(1.504N1f7, RoundUp) === 2N1f7
        @test round(1.504N1f7, RoundDown) === 1N1f7
        @test round(Int, 1.504N1f7, RoundNearest) === 2
        @test round(Int, 1.504N1f7, RoundToZero) === 1
        @test round(Int, 1.504N1f7, RoundUp) === 2
        @test round(Int, 1.504N1f7, RoundDown) === 1
    end
end

@testset "approx" begin
    @testset "approx $T" for T in FixedPointNumbers.UF
        xs = typemin(T):eps(T):typemax(T)-eps(T)
        @test all(x -> x ≈ x + eps(T), xs)
        @test all(x -> x + eps(T) ≈ x, xs)
        @test !any(x -> x - eps(T) ≈ x + eps(T), xs)
    end
end

@testset "low-level arithmetic" begin
    @test !(N0f8(0.5) < N0f8(0.5))
    @test N0f8(0.5) <= N0f8(0.5)

    @test div(reinterpret(N0f8, 0x10), reinterpret(N0f8, 0x02)) == fld(reinterpret(N0f8, 0x10), reinterpret(N0f8, 0x02)) == 8
    @test div(reinterpret(N0f8, 0x0f), reinterpret(N0f8, 0x02)) == fld(reinterpret(N0f8, 0x0f), reinterpret(N0f8, 0x02)) == 7
    @test Base.fld1(reinterpret(N0f8, 0x10), reinterpret(N0f8, 0x02)) == 8
    @test Base.fld1(reinterpret(N0f8, 0x0f), reinterpret(N0f8, 0x02)) == 8
    @test mod(reinterpret(N0f8, 0x10), reinterpret(N0f8, 0x02)) == rem(reinterpret(N0f8, 0x10), reinterpret(N0f8, 0x02)) == 0
    @test mod(reinterpret(N0f8, 0x0f), reinterpret(N0f8, 0x02)) == rem(reinterpret(N0f8, 0x0f), reinterpret(N0f8, 0x02)) == reinterpret(N0f8, 0x01)
    @test mod1(reinterpret(N0f8, 0x10), reinterpret(N0f8, 0x02)) == reinterpret(N0f8, 0x02)
    @test mod1(reinterpret(N0f8, 0x0f), reinterpret(N0f8, 0x02)) == reinterpret(N0f8, 0x01)
    @test bswap(N0f8(0.5)) === N0f8(0.5)
    @test bswap(N0f16(0.5)) === reinterpret(N0f16, 0x0080)
    @test minmax(N0f8(0.8), N0f8(0.2)) === (N0f8(0.2), N0f8(0.8))
end

@testset "clamp" begin
    @test clamp(0.5N0f8, 0.2N0f8, 0.8N0f8) === 0.5N0f8
    @test clamp(0.5N0f8, 0.6N0f8, 0.8N0f8) === 0.6N0f8
    @test clamp(0.5N0f8, 0.2N0f8, 0.4N0f8) === 0.4N0f8
    @test clamp(0.5,      0.2N0f8, 0.8N0f8) === 0.5
    @test clamp(0.5f0,    0.6N0f8, 0.8N0f8) === 0.6f0
    @test clamp(0.5N0f16, 0.2N0f8, 0.4N0f8) === 0.4N0f16
    @test clamp(0.6N0f8, -Inf, Inf) === 0.6
    @test clamp(0.5,     N0f8) === 0.5N0f8
    @test clamp(-1.0f0,  N0f8) === 0.0N0f8
    @test clamp(2.0N1f7, N0f8) === 1.0N0f8
end

@testset "sign-related functions" begin
    @test_throws Exception signed(N0f8)
    @test_throws Exception signed(1N0f8)
    @test_throws Exception unsigned(N0f8)
    @test_throws Exception unsigned(1N0f8)
    @test_throws ArgumentError copysign(1N0f8, 0x1)
    @test_throws ArgumentError copysign(1N0f8, -1)
    @test_throws ArgumentError flipsign(1N0f8, 0x1)
    @test_throws ArgumentError flipsign(1N0f8, -1)
    @test_throws ArgumentError sign(0N0f8)
    @test signbit(1N0f8) === false
end

@testset "unit range" begin
    @test length(N0f8(0):N0f8(1)) == 2
    @test length(N0f8(1):N0f8(0)) == 0
    @test isempty(N0f8(1):N0f8(0))
    @test collect(N0f8(0):N0f8(1)) == N0f8[0, 1]
    @test length(0.5N1f7:1.504N1f7) == 2
    @test length(N7f1(0):N7f1(255)) == 256
    NIntW = Normed{UInt,bitwidth(UInt)}
    @test length(NIntW(0):NIntW(1)) == 2
    NInt1 = Normed{UInt,1}
    @test length(NInt1(0):typemax(NInt1)-oneunit(NInt1)) == typemax(UInt)
    @test_throws OverflowError length(NInt1(0):typemax(NInt1))
    @test Base.unsafe_length(NInt1(0):typemax(NInt1)) == 0  # overflow
    N64f64 = Normed{UInt128,64}
    @test_broken length(N64f64(0):typemax(N64f64)) == UInt128(typemax(UInt64)) + 1
    @test length(N1f63(2):N1f63(0)) == 0
end

@testset "step range" begin
    counter = 0
    for x in N0f8(0):eps(N0f8):N0f8(1)
        counter += 1
    end
    @test counter == 256
    @test length(N0f8(0):eps(N0f8):N0f8(1)) == 256
    r = reinterpret(N0f8, 0x01):reinterpret(N0f8, 0x01):reinterpret(N0f8, UInt8(48))
    @test length(r) == 48
    NInt1 = Normed{UInt,1}
    @test length(NInt1(0):NInt1(1):typemax(NInt1)-oneunit(NInt1)) == typemax(UInt)
    @test_throws OverflowError length(NInt1(0):NInt1(1):typemax(NInt1))
end

@testset "predicates" begin
    @test isfinite(1N8f8)
    @test !isnan(1N8f8)
    @test !isinf(1N8f8)

    @testset "isinteger" begin
        for T in (UInt8, UInt16)
            @testset "isinteger(::Normed{$T,$f})" for f = 1:bitwidth(T)
                N = Normed{T,f}
                xs = typemin(N):eps(N):typemax(N)
                @test all(x -> isinteger(x) == isinteger(float(x)), xs)
            end
        end
        for T in (UInt32, UInt64)
            @testset "isinteger(::Normed{$T,$f})" for f = 1:bitwidth(T)
                N = Normed{T,f}
                if f == 1
                    @test isinteger(zero(N)) & isinteger(oneunit(N))
                else
                    @test !isinteger(oneunit(N) - eps(N)) & isinteger(oneunit(N))
                end
            end
        end
    end
end

@testset "Promotion within Normed" begin
    @test @inferred(promote(N0f8(0.2), N0f8(0.8))) ===
        (N0f8(0.2), N0f8(0.8))
    @test @inferred(promote(Normed{UInt16,3}(0.2), Normed{UInt8,3}(0.86))) ===
        (Normed{UInt16,3}(0.2), Normed{UInt16,3}(0.86))
    @test @inferred(promote(Normed{UInt8,7}(0.197), Normed{UInt8,4}(0.8))) ===
        (Normed{UInt16,7}(0.197), Normed{UInt16,7}(0.8))

    @test Normed{UInt16,16}(1)   == Normed{UInt8,8}(1)
    @test Normed{UInt16,16}(0.2) == Normed{UInt8,8}(0.2)
    @test Normed{UInt16,8}(1)   == Normed{UInt8,8}(1)
    @test Normed{UInt16,8}(0.2) == Normed{UInt8,8}(0.2)
    @test Normed{UInt16,16}(1)   == Normed{UInt8,6}(1)
    @test Normed{UInt16,16}(0.20635) == Normed{UInt8,6}(0.20635)
    @test Normed{UInt16,4}(1)   == Normed{UInt8,6}(1)
    @test Normed{UInt16,4}(0.2) == Normed{UInt8,6}(0.2)

    @test promote_type(N0f8,Float32,Int) == Float32
    @test promote_type(N0f8,Int,Float32) == Float32
    @test promote_type(Int,N0f8,Float32) == Float32
    @test promote_type(Int,Float32,N0f8) == Float32
    @test promote_type(Float32,Int,N0f8) == Float32
    @test promote_type(Float32,N0f8,Int) == Float32
    @test promote_type(N0f8,N1f7,N2f6,N3f5,N4f4,N5f3) == Normed{UInt128,8}
end

@testset "show" begin
    iob = IOBuffer()
    n0f8 = reinterpret(N0f8, 0xaa)
    show(iob, n0f8)
    str = String(take!(iob))
    @test str == "0.667N0f8"
    @test eval(Meta.parse(str)) === n0f8

    n16f16 = reinterpret(N16f16, 0xaaaaaaaa)
    show(iob, n16f16)
    str = String(take!(iob))
    @test str == "43691.33333N16f16"
    @test eval(Meta.parse(str)) === n16f16

    show(IOContext(iob, :compact=>true), n16f16)
    @test String(take!(iob)) == "43691.3"

    show(IOContext(iob, :compact=>true, :typeinfo=>N16f16), n16f16)
    @test String(take!(iob)) == "43691.3"

    show(IOContext(iob, :compact=>true, :typeinfo=>Normed), n16f16)
    @test String(take!(iob)) == "43691.3"

    show(IOContext(iob, :typeinfo=>N16f16), n16f16)
    @test String(take!(iob)) == "43691.33333N16f16" # TODO: Consider removing suffix (issue #188)

    show(IOContext(iob, :typeinfo=>Normed), n16f16)
    @test String(take!(iob)) == "43691.33333N16f16"

    show(iob, Normed{UInt128,64}(1.2345e6))
    @test_broken String(take!(iob)) == "Normed{UInt128,64}(1.2345e6)" # "N64f64" is not defined
end

@testset "summary" begin
    a = N0f8[0.2, 0.4]
    aa = Normed[0.2N0f8 0.4N0f16]

    if VERSION >= v"1.6.0-DEV.356"
        @test_broken summary(a) == "2-element Vector{N0f8}"
        @test_broken summary(view(a, 1:2)) == "2-element view(::Vector{N0f8}, 1:2) with eltype N0f8"
        @test_broken summary(aa) == "1×2 Matrix{Normed}"
    else
        @test summary(a) == "2-element Array{N0f8,1} with eltype Normed{UInt8,8}"
        @test summary(view(a, 1:2)) == "2-element view(::Array{N0f8,1}, 1:2) with eltype Normed{UInt8,8}"
        @test_broken summary(aa) == "1×2 Array{Normed,2}"
    end
end

@testset "scaledual" begin
    a = rand(UInt8, 10)
    af8 = reinterpret(N0f8, a)
    b = 0.5

    # LHSs of the following `@test`s with `af8` can be slightly more accurate
    bd, eld = scaledual(b, af8[1])
    @test b*af8[1] ≈ bd*eld rtol=1e-15
    bd, ad = scaledual(b, af8)
    @test b*af8 ≈ bd*ad rtol=1e-15

    bd, eld = scaledual(b, a[1])
    @test b*a[1] == bd*eld
    bd, ad = scaledual(b, a)
    @test b*a == bd*ad

    bd, eld = scaledual(Float64, af8[1])
    @test 1.0*af8[1] ≈ bd*eld rtol=1e-15
    bd, ad = scaledual(Float64, af8)
    @test 1.0*af8 ≈ bd*ad rtol=1e-15

    bd, eld = scaledual(Float64, a[1])
    @test 1.0*a[1] == bd*eld
    bd, ad = scaledual(Float64, a)
    @test 1.0*a == bd*ad
end

@testset "reductions" begin
    a = N0f8[reinterpret(N0f8, 0xff), reinterpret(N0f8, 0xff)]
    @test sum(a) == 2.0
    @test sum(a, dims=1) == [2.0]

    a = N2f14[3.2, 2.4]
    acmp = Float64(a[1])*Float64(a[2])
    @test prod(a) == acmp
    @test prod(a, dims=1) == [acmp]
end

@testset "reductions, Statistics" begin
    a = N0f8[reinterpret(N0f8, 0x80), reinterpret(N0f8, 0x40)]
    af = FixedPointNumbers.Treduce.(a)
    @test mean(a) === mean(af)
    @test std(a)  === std(af)
    @test var(a)  === var(af)
    m = mean(a)
    @test stdm(a, m) === stdm(af, m)
    @test varm(a, m) === varm(af, m)
end

@testset "rand" begin
    for T in (Normed{UInt8,8}, Normed{UInt8,6},
              Normed{UInt16,16}, Normed{UInt16,14},
              Normed{UInt32,32}, Normed{UInt32,30},
              Normed{UInt64,64}, Normed{UInt64,62})
        a = rand(T)
        @test isa(a, T)
        a = rand(T, (3, 5))
        @test ndims(a) == 2 && eltype(a) == T
        @test size(a) == (3,5)
    end
end

@testset "Overflow with Float16" begin
    @test N0f16(Float16(1.0)) === N0f16(1.0)
    @test Float16(1.0) % N0f16 === N0f16(1.0)
end

@testset "disambiguation constructors" begin
    @test_throws ArgumentError Normed{UInt32,16}('a')
    @test_throws InexactError  Normed{UInt32,16}(complex(1.0, 1.0))
    @test Normed{UInt32,16}(complex(1.0, 0.0))        == 1
    @test Normed{UInt32,16}(Base.TwicePrecision(1.0, 0.0)) == 1
end
