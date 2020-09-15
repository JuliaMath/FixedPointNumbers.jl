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

@testset "construction using type suffix" begin
    @test 0.635N0f8  ===  N0f8(0.635)
    @test 0.635N6f10 === N6f10(0.635)
    @test 0.635N4f12 === N4f12(0.635)
    @test 0.635N2f14 === N2f14(0.635)
    @test 0.635N0f16 === N0f16(0.635)

    @test_throws ArgumentError 1.1N0f8

    @test wrapping_mul(0.635, N0f8) === N0f8(0.635)
    @test wrapping_mul(1.635, N0f8) === N0f8((1.635 * 255 - 256) / 255)
    @test saturating_mul(0.635, N0f8) === N0f8(0.635)
    @test saturating_mul(1.635, N0f8) === N0f8(1.0)
    @test checked_mul(0.635, N0f8) === N0f8(0.635)
    @test_throws ArgumentError checked_mul(1.635, N0f8)
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
    @test occursin("Normed{UInt128,$(SP)100} is a 128-bit type representing 2^128 values", msg)
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

    test_floattype(Normed)
end

@testset "conversions from float" begin
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

    test_convert_from_nan(Normed)
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
            float_err = zero(Tf)
            for i = typemin(T):typemax(T)
                f_expected = Tf(i / BigFloat(FixedPointNumbers.rawone(N)))
                isinf(f_expected) && break # for Float16(::Normed{UInt16,1})
                f_actual = Tf(reinterpret(N, i))
                float_err += abs(f_actual - f_expected)
            end
            @test float_err == 0
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
    test_rem_type(Normed)
    test_rem_nan(Normed)

    @test  N0f8(0.2) % N0f8  === N0f8(0.2)
    @test N2f14(1.2) % N0f16 === N0f16(0.20002)
    @test N2f14(1.2) % N0f8  === N0f8(0.196)

    @test ( 1.5 % N0f8).i == round(Int,  1.5*255) % UInt8
    @test (-0.3 % N0f8).i == round(Int, -0.3*255) % UInt8

    @test (65.2 % N6f10).i == round(Int, 65.2*1023) % UInt16
    @test (-0.3 % N6f10).i == round(Int, -0.3*1023) % UInt16

    @test 1 % N0f8 === N0f8(1)
    @test 2 % N0f8 === N0f8(0.996)

    # issue #150
    @test all(f -> 1.0f0 % Normed{UInt32,f} == oneunit(Normed{UInt32,f}), 1:32)
    @test all(f -> 1.0e0 % Normed{UInt64,f} == oneunit(Normed{UInt64,f}), 1:64)

    # issue #211
    @test big"1.2" % N0f8 === 0.196N0f8
    @test reinterpret(BigFloat(0x0_01234567_89abcdef) % N63f1) === 0x01234567_89abcdef
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
    test_neg(Normed)
end

@testset "abs" begin
    for N in target(Normed; ex = :thin)
        @test   wrapping_abs(typemax(N)) === typemax(N)
        @test saturating_abs(typemax(N)) === typemax(N)
        @test    checked_abs(typemax(N)) === typemax(N)

        @test   wrapping_abs(typemin(N)) === typemin(N)
        @test saturating_abs(typemin(N)) === typemin(N)
        @test    checked_abs(typemin(N)) === typemin(N)
    end
    test_abs(Normed)
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
    test_add(Normed)
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
    test_sub(Normed)
end

@testset "mul" begin
    for N in target(Normed; ex = :thin)
        @test   wrapping_mul(typemax(N), zero(N)) === zero(N)
        @test saturating_mul(typemax(N), zero(N)) === zero(N)
        @test    checked_mul(typemax(N), zero(N)) === zero(N)

        @test   wrapping_mul(one(N), typemax(N)) === typemax(N)
        @test saturating_mul(one(N), typemax(N)) === typemax(N)
        @test    checked_mul(one(N), typemax(N)) === typemax(N)

        @test   wrapping_mul(typemax(N), typemax(N)) === big(typemax(N))^2 % N
        @test saturating_mul(typemax(N), typemax(N)) === typemax(N)
        if typemax(N) != 1
            @test_throws OverflowError checked_mul(typemax(N), typemax(N))
        end
    end
    test_mul(Normed)
end

@testset "fdiv" begin
    for N in target(Normed; ex = :thin)
        @test   wrapping_fdiv(typemax(N), typemax(N)) === one(N)
        @test saturating_fdiv(typemax(N), typemax(N)) === one(N)
        @test    checked_fdiv(typemax(N), typemax(N)) === one(N)

        @test   wrapping_fdiv(zero(N), eps(N)) === zero(N)
        @test saturating_fdiv(zero(N), eps(N)) === zero(N)
        @test    checked_fdiv(zero(N), eps(N)) === zero(N)

        @test   wrapping_fdiv(typemax(N), eps(N)) === (floattype(N))(typemax(rawtype(N))) % N
        @test saturating_fdiv(typemax(N), eps(N)) === typemax(N)
        @test_throws OverflowError checked_fdiv(typemax(N), eps(N))

        @test   wrapping_fdiv(zero(N), zero(N)) === zero(N)
        @test saturating_fdiv(zero(N), zero(N)) === zero(N)
        @test_throws DivideError checked_fdiv(zero(N), zero(N))

        @test   wrapping_fdiv(eps(N), zero(N)) === zero(N)
        @test saturating_fdiv(eps(N), zero(N)) === typemax(N)
        @test_throws DivideError checked_fdiv(eps(N), zero(N))
    end
    test_fdiv(Normed)
end

@testset "div/cld/fld" begin
    for N in target(Normed; ex = :thin)
        nm, nz, ne = typemax(N), zero(N), eps(N)
        T = rawtype(N)
        @test   wrapping_div(nm, nm) ===   wrapping_fld(nm, nm) ===   wrapping_cld(nm, nm) === one(T)
        @test saturating_div(nm, nm) === saturating_fld(nm, nm) === saturating_cld(nm, nm) === one(T)
        @test    checked_div(nm, nm) ===    checked_fld(nm, nm) ===    checked_cld(nm, nm) === one(T)

        @test   wrapping_div(nz, ne) ===   wrapping_fld(nz, ne) ===   wrapping_cld(nz, ne) === zero(T)
        @test saturating_div(nz, ne) === saturating_fld(nz, ne) === saturating_cld(nz, ne) === zero(T)
        @test    checked_div(nz, ne) ===    checked_fld(nz, ne) ===    checked_cld(nz, ne) === zero(T)

        @test   wrapping_div(nm, ne) ===   wrapping_fld(nm, ne) ===   wrapping_cld(nm, ne) === typemax(T)
        @test saturating_div(nm, ne) === saturating_fld(nm, ne) === saturating_cld(nm, ne) === typemax(T)
        @test    checked_div(nm, ne) ===    checked_fld(nm, ne) ===    checked_cld(nm, ne) === typemax(T)

        @test   wrapping_div(nz, nz) ===   wrapping_fld(nz, nz) ===   wrapping_cld(nz, nz) === zero(T)
        @test saturating_div(nz, nz) === saturating_fld(nz, nz) === saturating_cld(nz, nz) === zero(T)
        @test_throws DivideError checked_div(nz, nz)
        @test_throws DivideError checked_fld(nz, nz)
        @test_throws DivideError checked_cld(nz, nz)

        @test   wrapping_div(ne, nz) ===   wrapping_fld(ne, nz) ===   wrapping_cld(ne, nz) === zero(T)
        @test saturating_div(ne, nz) === saturating_fld(ne, nz) === saturating_cld(ne, nz) === typemax(T)
        @test_throws DivideError checked_div(ne, nz)
        @test_throws DivideError checked_fld(ne, nz)
        @test_throws DivideError checked_cld(ne, nz)

        @test wrapping_div(ne, nm) === saturating_div(ne, nm) === checked_div(ne, nm) === zero(T)
        @test wrapping_fld(ne, nm) === saturating_fld(ne, nm) === checked_fld(ne, nm) === zero(T)
        @test wrapping_cld(ne, nm) === saturating_cld(ne, nm) === checked_cld(ne, nm) === one(T)
    end
    test_div(Normed)
    test_div_3arg(Normed)
end

@testset "rem/mod" begin
    for N in target(Normed; ex = :thin)
        nm, nz, ne = typemax(N), zero(N), eps(N)
        T = rawtype(N)
        @test   wrapping_rem(nm, nm) ===   wrapping_mod(nm, nm) === nz
        @test saturating_rem(nm, nm) === saturating_mod(nm, nm) === nz
        @test    checked_rem(nm, nm) ===    checked_mod(nm, nm) === nz

        @test   wrapping_rem(nz, ne) ===   wrapping_mod(nz, ne) === nz
        @test saturating_rem(nz, ne) === saturating_mod(nz, ne) === nz
        @test    checked_rem(nz, ne) ===    checked_mod(nz, ne) === nz

        @test   wrapping_rem(nm, ne) ===   wrapping_mod(nm, ne) === nz
        @test saturating_rem(nm, ne) === saturating_mod(nm, ne) === nz
        @test    checked_rem(nm, ne) ===    checked_mod(nm, ne) === nz

        @test   wrapping_rem(nz, nz) ===   wrapping_mod(nz, nz) === nz
        @test saturating_rem(nz, nz) === saturating_mod(nz, nz) === nz
        @test_throws DivideError checked_rem(nz, nz)
        @test_throws DivideError checked_mod(nz, nz)

        @test   wrapping_rem(ne, nz) ===   wrapping_mod(ne, nz) === ne
        @test saturating_rem(ne, nz) === saturating_mod(ne, nz) === ne
        @test_throws DivideError checked_rem(ne, nz)
        @test_throws DivideError checked_mod(ne, nz)

        @test wrapping_rem(ne, nm) === saturating_rem(ne, nm) === checked_rem(ne, nm) === ne
        @test wrapping_mod(ne, nm) === saturating_mod(ne, nm) === checked_mod(ne, nm) === ne
    end
    test_rem(Normed)
    test_rem_3arg(Normed)

    @test_throws OverflowError rem(0.5N0f8, 1N0f8, RoundUp)
    @test saturating_rem(0.5N0f8, 1N0f8, RoundUp) === zero(N0f8)
end

@testset "fld1/mod1" begin
    test_fld1_mod1(Normed)
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
    test_isapprox(Normed)

    @test isapprox(typemin(N0f8), typemax(N0f8), rtol=1.0)
    @test !isapprox(zero(N0f8), typemax(N0f8), rtol=0.9)
    @test isapprox(zero(N0f8), eps(N0f8), rtol=1e-6) # atol = eps(N0f8)
    @test !isapprox(eps(N0f8), zero(N0f8), rtol=1e-6, atol=1e-6)
    @test !isapprox(0.66N6f2, 1.0N6f2, rtol=0.3, atol=0) # 1.0 * 0.3 < eps(N6f2)

    @test isapprox(eps(N8f8), eps(N0f8), rtol=1e-6)
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

    test_clamp_nan(Normed)
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
        test_isinteger(Normed)
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
    test_rand(Normed)
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
    @test String(take!(iob)) == "Normed{UInt128,$(SP)64}(1.2345e6)"
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
