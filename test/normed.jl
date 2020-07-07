include("common.jl")

function symbol_to_inttype(::Type{Normed}, s::Symbol)
    d = Dict(:i8 => UInt8, :i16 => UInt16, :i32 => UInt32, :i64 => UInt64, :i128 => UInt128)
    d[s]
end

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

@testset "limits and identities" begin
    @testset "$N" for N in target(Normed)
        T, f = rawtype(N), nbitsfrac(N)
        @test zero(N) == 0
        @test one(N) == 1
        @test one(N) * oneunit(N) == oneunit(N)
        @test typemin(N) == 0
        @test typemax(N) == typemax(T)//(big"2"^f - 1)
        @test floatmin(N) === eps(N) == 1//(big"2"^f - 1)
        @test floatmax(N) === typemax(N)
        @test eps(zero(N)) === eps(typemax(N))
        @test sizeof(N) == sizeof(T)
    end
end

@testset "inexactness" begin
    # TODO: change back to InexactError when it allows message strings
    @test_throws ArgumentError N0f8(2)
    @test_throws ArgumentError N0f8(0xff)
    @test_throws ArgumentError N0f16(2)
    @test_throws ArgumentError N0f16(0xff)
    @test_throws ArgumentError N0f16(0xffff)
    @test_throws ArgumentError convert(N0f8,  typemax(N6f10))
    @test_throws ArgumentError convert(N0f16, typemax(N6f10))
    @test_throws ArgumentError convert(Normed{UInt128,100}, 10^9)

    ret = @test_throws ArgumentError N0f8(255)
    msg = ret.value.msg
    @test occursin("N0f8 is an 8-bit type representing 256 values from 0.0 to 1.0;", msg)
    ret = @test_throws ArgumentError convert(Normed{UInt128,100}, 10.0^9)
    msg = ret.value.msg
    @test occursin("Normed{UInt128,100} is a 128-bit type representing 2^128 values", msg)
end

@testset "disambiguation constructors" begin
    @test_throws ArgumentError Normed{UInt32,16}('a')
    @test_throws InexactError  Normed{UInt32,16}(complex(1.0, 1.0))
    @test Normed{UInt32,16}(complex(1.0, 0.0)) == 1
    @test Normed{UInt32,16}(Base.TwicePrecision(1.0, 0.0)) == 1
end

@testset "conversion" begin
    x = N0f8(0.5)
    @test convert(N0f8, x) === x

    @test convert(N0f8,  1.1/typemax(UInt8)) === eps(N0f8)
    @test convert(N0f8,  1.1f0/typemax(UInt8)) === eps(N0f8)
    @test convert(N6f10, 1.1/typemax(UInt16)*64) === eps(N6f10)
    @test convert(N4f12, 1.1/typemax(UInt16)*16) === eps(N4f12)
    @test convert(N2f14, 1.1/typemax(UInt16)*4)  === eps(N2f14)
    @test convert(N0f16, 1.1/typemax(UInt16))    === eps(N0f16)
    @test convert(N16f16, 1.1/typemax(UInt32)*2^16) === eps(N16f16)
    @test convert(N61f3,  1.1/typemax(UInt64)*UInt64(2)^61) === eps(N61f3)
    @test convert(Normed{UInt128,7}, 1.1/typemax(UInt128)*UInt128(2)^121) === eps(Normed{UInt128,7})

    @test convert(N0f8, Base.TwicePrecision(1.0)) === 1N0f8

    @test convert(N0f16, one(N0f8)) === one(N0f16)
    @test convert(N0f16, N0f8(0.5)) === reinterpret(N0f16, 0x8080)
    @test convert(N9f7, N1f7(0.504)) === N9f7(0.504)

    # avoiding overflow with Float16
    @test N0f16(Float16(1.0)) === N0f16(1.0)
    @test Float16(1.0) % N0f16 === N0f16(1.0)
end

@testset "bool conversions" begin
    @testset "$N to/from Bool" for N in target(Normed)
        @test convert(Bool, zero(N)) === false
        @test convert(Bool, oneunit(N))  === true
        eps(N) < 1 && @test_throws InexactError convert(Bool, convert(N, 0.2))
        @test convert(N, true) === oneunit(N)
        @test convert(N, false) === zero(N)
    end
    @test Bool(1N0f8) === true
end

@testset "integer conversions" begin
    @testset "$N to/from integer" for N in target(Normed)
        @test convert(Int, oneunit(N)) === 1
        @test convert(Integer, oneunit(N)) === oneunit(rawtype(N))
        @test convert(N, 1) === oneunit(N)
        @test convert(N, 0x0) === zero(N)
    end
    @test convert(UInt, 1N1f7) === UInt(1)
    @test_throws InexactError convert(Integer, 0.5N1f7)
    @test_throws InexactError convert(Int8, 256N8f8)
end

@testset "rational conversions" begin
    @testset "$N to/from rational" for N in target(Normed)
        @test convert(Rational, oneunit(N)) == 1//1
        @test convert(Rational{Int}, zero(N)) === 0//1
        @test convert(N, 1//1) === oneunit(N)
    end
    @test convert(N0f8, 1//255) === eps(N0f8)
    @test convert(N0f8, Rational{Int8}(3//5)) === N0f8(3/5)
    @test convert(N0f8, Rational{UInt8}(3//5)) === N0f8(3/5)
    @test_throws ArgumentError convert(N0f8, typemax(Rational{UInt8}))

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
    @test convert(BigFloat, eps(N0f8))::BigFloat == 1 / big"255"

    @test big(N7f1) === BigFloat # !== BigInt
    @test big(0.5N4f4)::BigFloat == 8 / big"15"
end

@testset "float/floattype" begin
    @test float(0.8N4f4) === 0.8f0
    @test float(0.8N20f12) === 0.8
    @test float(0.8N8f24) === 0.8
    @test float(1N11f53)::BigFloat == big"1.0"

    @testset "floattype($N)" for N in target(Normed, :i8, :i16, :i32, :i64; ex = :heavy)
        @test typemax(N) <= maxintfloat(floattype(N))
    end
end

@testset "conversion from float" begin
    # issue 102
    for Tf in (Float16, Float32, Float64)
        @testset "$N(::$Tf)" for N in target(Normed)
            T, f = rawtype(N), nbitsfrac(N)
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
        @testset "$Tf(::$N)" for N in target(Normed, :i8, :i16)
            T = rawtype(N)
            float_err = 0.0
            for i = typemin(T):typemax(T)
                f_expected = Tf(i / BigFloat(FixedPointNumbers.rawone(N)))
                isinf(f_expected) && break # for Float16(::Normed{UInt16,1})
                f_actual = Tf(reinterpret(N, i))
                float_err += abs(f_actual - f_expected)
            end
            @test float_err == 0.0
        end
        @testset "$Tf(::$N)" for N in target(Normed, :i32, :i64, :i128)
            T = rawtype(N)
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

@testset "type modulus" begin
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

@testset "arithmetic" begin
    @testset "$N arithmetic" for N in target(Normed; ex = :light)
        x = N(0x10,0)
        y = N(0x25,0)
        fx = float(x)
        fy = float(y)
        @test y > x
        @test y != x
        @test x+y === N(0x35,0)
        @test ((x+y)-y) === x
        @test ((x-y)+y) === x # wraparound
        fx*fy <= typemax(N) && @test (x*y)::N ≈ convert(N, fx*fy)
        @test (x/y)::N ≈ convert(N, fx/fy)
        fx^2 <= typemax(N) && @test (x^2)::N ≈ convert(N, fx^2)
        @test (x^2.1f0) ≈ fx^2.1f0
        @test (x^2.1) ≈ convert(Float64, x)^2.1
    end
end

@testset "neg" begin
    for N in target(Normed; ex = :thin)
        @test   wrapping_neg(typemin(N)) === zero(N)
        @test saturating_neg(typemin(N)) === zero(N)
        @test    checked_neg(typemin(N)) === zero(N)

        @test   wrapping_neg(typemax(N)) === eps(N)
        @test saturating_neg(typemax(N)) === zero(N)
        @test_throws OverflowError checked_neg(typemax(N))

        @test   wrapping_neg(eps(N)) === typemax(N)
        @test saturating_neg(eps(N)) === zero(N)
        @test_throws OverflowError checked_neg(eps(N))
    end
    for N in target(Normed, :i8; ex = :thin)
        xs = typemin(N):eps(N):typemax(N)
        fneg(x) = -float(x)
        @test all(x -> wrapping_neg(wrapping_neg(x)) === x, xs)
        @test all(x -> saturating_neg(x) === clamp(fneg(x), N), xs)
        @test all(x -> !(typemin(N) < fneg(x) < typemax(N)) ||
                       wrapping_neg(x) === checked_neg(x) === fneg(x) % N, xs)
    end
end

@testset "add" begin
    for N in target(Normed; ex = :thin)
        @test   wrapping_add(typemin(N), typemin(N)) === zero(N)
        @test saturating_add(typemin(N), typemin(N)) === zero(N)
        @test    checked_add(typemin(N), typemin(N)) === zero(N)

        @test   wrapping_add(typemax(N), eps(N)) ===   wrapping_add(eps(N), typemax(N)) === zero(N)
        @test saturating_add(typemax(N), eps(N)) === saturating_add(eps(N), typemax(N)) === typemax(N)
        @test_throws OverflowError checked_add(typemax(N), eps(N))
        @test_throws OverflowError checked_add(eps(N), typemax(N))

        @test   wrapping_add(zero(N), eps(N)) ===   wrapping_add(eps(N), zero(N)) === eps(N)
        @test saturating_add(zero(N), eps(N)) === saturating_add(eps(N), zero(N)) === eps(N)
        @test    checked_add(zero(N), eps(N)) ===    checked_add(eps(N), zero(N)) === eps(N)
    end
    for N in target(Normed, :i8; ex = :thin)
        xs = typemin(N):eps(N):typemax(N)
        xys = ((x, y) for x in xs, y in xs)
        fadd(x, y) = float(x) + float(y)
        @test all(((x, y),) -> wrapping_sub(wrapping_add(x, y), y) === x, xys)
        @test all(((x, y),) -> saturating_add(x, y) === clamp(fadd(x, y), N), xys)
        @test all(((x, y),) -> !(typemin(N) < fadd(x, y) < typemax(N)) ||
                               wrapping_add(x, y) === checked_add(x, y) === fadd(x, y) % N, xys)
    end
end

@testset "sub" begin
    for N in target(Normed; ex = :thin)
        @test   wrapping_sub(typemin(N), typemin(N)) === zero(N)
        @test saturating_sub(typemin(N), typemin(N)) === zero(N)
        @test    checked_sub(typemin(N), typemin(N)) === zero(N)

        @test   wrapping_sub(typemin(N), eps(N)) === typemax(N)
        @test saturating_sub(typemin(N), eps(N)) === typemin(N)
        @test_throws OverflowError checked_sub(typemin(N), eps(N))

        @test   wrapping_sub(eps(N), zero(N)) === eps(N)
        @test saturating_sub(eps(N), zero(N)) === eps(N)
        @test    checked_sub(eps(N), zero(N)) === eps(N)
    end
    for N in target(Normed, :i8; ex = :thin)
        xs = typemin(N):eps(N):typemax(N)
        xys = ((x, y) for x in xs, y in xs)
        fsub(x, y) = float(x) - float(y)
        @test all(((x, y),) -> wrapping_add(wrapping_sub(x, y), y) === x, xys)
        @test all(((x, y),) -> saturating_sub(x, y) === clamp(fsub(x, y), N), xys)
        @test all(((x, y),) -> !(typemin(N) < fsub(x, y) < typemax(N)) ||
                               wrapping_sub(x, y) === checked_sub(x, y) === fsub(x, y) % N, xys)
    end
end

@testset "div/fld1" begin
    @test div(reinterpret(N0f8, 0x10), reinterpret(N0f8, 0x02)) == fld(reinterpret(N0f8, 0x10), reinterpret(N0f8, 0x02)) == 8
    @test div(reinterpret(N0f8, 0x0f), reinterpret(N0f8, 0x02)) == fld(reinterpret(N0f8, 0x0f), reinterpret(N0f8, 0x02)) == 7
    @test fld1(reinterpret(N0f8, 0x10), reinterpret(N0f8, 0x02)) == 8
    @test fld1(reinterpret(N0f8, 0x0f), reinterpret(N0f8, 0x02)) == 8
end

@testset "rem/mod" begin
    @test mod(reinterpret(N0f8, 0x10), reinterpret(N0f8, 0x02)) == rem(reinterpret(N0f8, 0x10), reinterpret(N0f8, 0x02)) == 0
    @test mod(reinterpret(N0f8, 0x0f), reinterpret(N0f8, 0x02)) == rem(reinterpret(N0f8, 0x0f), reinterpret(N0f8, 0x02)) == reinterpret(N0f8, 0x01)
    @test mod1(reinterpret(N0f8, 0x10), reinterpret(N0f8, 0x02)) == reinterpret(N0f8, 0x02)
    @test mod1(reinterpret(N0f8, 0x0f), reinterpret(N0f8, 0x02)) == reinterpret(N0f8, 0x01)
end

@testset "rounding" begin
    for sym in (:i8, :i16, :i32, :i64)
        T = symbol_to_inttype(Normed, sym)
        rs = vcat([ oneunit(T) << b - oneunit(T) << 1 for b = 1:bitwidth(T)],
                  [ oneunit(T) << b - oneunit(T)      for b = 1:bitwidth(T)],
                  [ oneunit(T) << b                   for b = 2:bitwidth(T)-1])
        @testset "rounding $N" for N in target(Normed, sym)
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
    @testset "approx $N" for N in target(Normed, :i8, :i16; ex = :light)
        xs = typemin(N):eps(N):typemax(N)-eps(N)
        @test all(x -> x ≈ x + eps(N), xs)
        @test all(x -> x + eps(N) ≈ x, xs)
        @test !any(x -> x - eps(N) ≈ x + eps(N), xs)
    end
end

@testset "comparison" begin
    @test !(N0f8(0.5) < N0f8(0.5))
    @test N0f8(0.5) <= N0f8(0.5)

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

@testset "bitwise" begin
    x = N0f8(0b01010001, 0)
    @test ~x == N0f8(0b10101110, 0)
    @test -x == reinterpret(N0f8, 0xaf)

    @test bswap(N0f8(0.5)) === N0f8(0.5)
    @test bswap(N0f16(0.5)) === reinterpret(N0f16, 0x0080)
end

@testset "predicates" begin
    @test isfinite(1N8f8)
    @test !isnan(1N8f8)
    @test !isinf(1N8f8)

    @testset "isinteger" begin
        @testset "isinteger(::$N)" for N in target(Normed, :i8, :i16)
            xs = typemin(N):eps(N):typemax(N)
            @test all(x -> isinteger(x) == isinteger(float(x)), xs)
        end
        @testset "isinteger(::$N)" for N in target(Normed, :i32, :i64, :i128)
            if nbitsfrac(N) == 1
                @test isinteger(zero(N)) & isinteger(oneunit(N))
            else
                @test !isinteger(oneunit(N) - eps(N)) & isinteger(oneunit(N))
            end
        end
    end
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
    @testset "rand(::$N)" for N in target(Normed; ex = :thin)
        @test isa(rand(N), N)
        a = rand(N, (3, 5))
        @test ndims(a) == 2 && eltype(a) === N
        @test size(a) == (3,5)
    end
    @test rand(MersenneTwister(1234), N0f8) === 0.925N0f8
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

    @test @inferred(promote_type(N0f8, Float64)) === Float64
    @test @inferred(promote_type(Float32, N8f24)) === Float64

    @test @inferred(promote_type(N0f8, Int8)) === Float32
    @test @inferred(promote_type(Int128, N8f24)) === Float64

    @test @inferred(promote_type(N0f16, Rational{Int8})) === Rational{Int8}

    @test @inferred(promote_type(N0f8, Float32, Int)) === Float32
    @test @inferred(promote_type(N0f8, Int, Float32)) === Float32
    @test @inferred(promote_type(Int, N0f8, Float32)) === Float32
    @test @inferred(promote_type(Int, Float32, N0f8)) === Float32
    @test @inferred(promote_type(Float32, Int, N0f8)) === Float32
    @test @inferred(promote_type(Float32, N0f8, Int)) === Float32

    @test @inferred(promote_type(N0f8,N1f7,N2f6,N3f5,N4f4,N5f3)) === Normed{UInt128,8}

    @test @inferred(promote_type(N0f8, Q0f31)) === Float64
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
    @test String(take!(iob)) == "43691.33333"

    show(IOContext(iob, :typeinfo=>Normed), n16f16)
    @test String(take!(iob)) == "43691.33333N16f16"

    show(iob, Normed{UInt128,64}(1.2345e6))
    @test String(take!(iob)) == "Normed{UInt128,64}(1.2345e6)"
end

@testset "summary" begin
    a = N0f8[0.2, 0.4]
    aa = Normed[0.2N0f8 0.4N0f16]

    if VERSION >= v"1.6.0-DEV.356"
        @test summary(a) == "2-element Vector{N0f8}"
        @test summary(view(a, 1:2)) == "2-element view(::Vector{N0f8}, 1:2) with eltype N0f8"
        @test summary(aa) == "1×2 Matrix{Normed}"
    else
        @test summary(a) == "2-element Array{N0f8,1} with eltype Normed{UInt8,8}"
        @test summary(view(a, 1:2)) == "2-element view(::Array{N0f8,1}, 1:2) with eltype Normed{UInt8,8}"
        @test summary(aa) == "1×2 Array{Normed,2}"
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
