include("common.jl")

function symbol_to_inttype(::Type{Fixed}, s::Symbol)
    d = Dict(:i8 => Int8, :i16 => Int16, :i32 => Int32, :i64 => Int64, :i128 => Int128)
    d[s]
end

function test_op(fun::Fun, fx::F, fy::F, fxf, fyf, tol) where {Fun, F}
    # Make sure that the result is representable
    zf = fun(fxf, fyf)
    typemin(F) <= zf <= typemax(F) || return nothing
    z = fun(fx, fy)
    @assert abs(z - convert(F, zf)) <= tol
    @assert abs(convert(Float64, z) - zf) <= tol
end

function test_fixed(::Type{F}) where {F}
    tol = Float64(eps(F))
    v = [-10:0.01:10; -180:.01:-160; 160:.01:180]
    # Ignore values outside the representable range
    values = filter(x -> typemin(F) < x <= typemax(F), v)
    for x in values
        fx = convert(F, x)
        fxf = convert(Float64, fx)

        @test convert(F, convert(Float64, fx)) === fx
        fx != typemin(F) && @test convert(F, convert(Float64, -fx)) === -fx
        fx != typemin(F) && @test convert(Float64, -fx) == -convert(Float64, fx)

        rx = convert(Rational{BigInt}, fx)
        @assert isequal(fx, rx) == isequal(hash(fx), hash(rx))

        for y in values
            fy = convert(F, y)
            fyf = convert(Float64, fy)

            @assert fx==fy || x!=y
            @assert fx<fy  || (x + tol)>=y
            @assert fx<=fy || x>y

            test_op(+, fx, fy, fxf, fyf, tol)
            test_op(-, fx, fy, fxf, fyf, tol)
            test_op(*, fx, fy, fxf, fyf, tol)
            fy != 0 && test_op(/, fx, fy, fxf, fyf, tol)

            @assert isequal(fx, fy) === isequal(hash(fx), hash(fy))
        end
    end
end

@testset "test_fixed" begin
    for F in target(Fixed, :i8, :i16, :i32; ex = :thin)
        test_fixed(F)
    end
end

@testset "domain of f" begin
    # TODO: change the upper limit
    @test_logs (:warn, r"`f=8` with raw type `T=Int8` will be removed") zero(Fixed{Int8,8})
    @test_throws DomainError zero(Fixed{Int8,-1})
    # @test_throws DomainError zero(Fixed{Int8,8})
    @test_throws DomainError zero(Fixed{Int8,9})
    # @test_throws DomainError zero(Fixed{Int16,16})
    @test_throws DomainError zero(Fixed{Int16,17})
end

@testset "reinterpret/bitstring" begin
    @test reinterpret(Q0f7, signed(0xa2)) === -0.734375Q0f7
    @test reinterpret(Q5f10, signed(0x00a2)) === 0.158203125Q5f10

    @test reinterpret(reinterpret(Q0f7, signed(0xa2))) === signed(0xa2)
    @test reinterpret(reinterpret(Q5f10, signed(0x00a2))) === signed(0x00a2)

    @test reinterpret(Int8, 0.5Q0f7) === signed(0x40)

    @test bitstring(reinterpret(Q0f7, signed(0xa2)))    === "10100010"
    @test bitstring(reinterpret(Q5f10, signed(0x00a2))) === "0000000010100010"
end

@testset "masks" begin
    @test FixedPointNumbers.intmask(0Q7f0) === signed(0xFF)
    @test FixedPointNumbers.intmask(0Q6f1) === signed(0xFE)
    @test FixedPointNumbers.intmask(0Q1f6) === signed(0xC0)
    @test FixedPointNumbers.intmask(0Q0f7) === signed(0x80)

    @test FixedPointNumbers.fracmask(0Q7f0) === signed(0x00)
    @test FixedPointNumbers.fracmask(0Q6f1) === signed(0x01)
    @test FixedPointNumbers.fracmask(0Q1f6) === signed(0x3F)
    @test FixedPointNumbers.fracmask(0Q0f7) === signed(0x7F)
end

@testset "limits and identities" begin
    @testset "$F" for F in target(Fixed)
        T, f = rawtype(F), nbitsfrac(F)
        @test zero(F) == 0
        f < bitwidth(T) - 1 && @test one(F) == 1
        f < bitwidth(T) - 1 && @test one(F) * oneunit(F) == oneunit(F)
        @test typemin(F) == typemin(T) >> f
        if T === Int128
            @test typemax(F) * big"2.0"^f == typemax(T) # force promotion to BigFloat due to lack of PR #207
        else
            @test typemax(F) == typemax(T)//big"2"^f
        end
        @test floatmin(F) === eps(F) == 2.0^-f # issue #79
        @test floatmax(F) === typemax(F)
        @test eps(zero(F)) === eps(typemax(F))
        @test sizeof(F) == sizeof(T)
    end
end

@testset "inexactness" begin
    # TODO: change back to InexactError when it allows message strings
    @test_throws ArgumentError Q0f7(-2)
    @test_throws ArgumentError one(Q0f15)
    @test_throws ArgumentError oneunit(Q0f31)
    @test_throws ArgumentError one(Fixed{Int8,8}) # TODO: remove this at end of its support

    @test_throws ArgumentError convert(Q0f7, 0.999)
    @test_throws ArgumentError convert(Q0f7, 1.0)
    @test_throws ArgumentError convert(Q0f7, 1)
    @test_throws ArgumentError convert(Q0f7, 2)

    ret = @test_throws ArgumentError Q0f7(127)
    msg = ret.value.msg
    @test occursin("Q0f7 is an 8-bit type representing 256 values from -1.0 to 0.992;", msg)
    ret = @test_throws ArgumentError convert(Fixed{Int128,100}, 10.0^9)
    msg = ret.value.msg
    @test occursin("Fixed{Int128,$(SP)100} is a 128-bit type representing 2^128 values", msg)
end

@testset "disambiguation constructors" begin
    @test_throws ArgumentError Fixed{Int32,16}('a')
    @test_throws InexactError  Fixed{Int32,16}(complex(1.0, 1.0))
    @test Fixed{Int32,16}(complex(1.0, 0.0)) == 1
    @test Fixed{Int32,16}(Base.TwicePrecision(1.0, 0.0)) == 1
end

@testset "conversion" begin
    @test isapprox(convert(Fixed{Int8,7}, 0.8), 0.797, atol=0.001)
    @test isapprox(convert(Fixed{Int8,7}, 0.9), 0.898, atol=0.001)

    @test convert(Q0f7, -128.5/128) == -1
    @test convert(Q0f7, -0.75f0) == -0.75
    @test convert(Q0f7, Float16(-0.75)) == -0.75
    @test convert(Q0f7, BigFloat(-0.75)) == -0.75
    @test_throws ArgumentError convert(Q0f7, BigFloat(127.5/128))

    @test convert(Q2f5, -1//2) === -0.5Q2f5
    @test convert(Q1f6, Rational{Int8}(-3//4)) === -0.75Q1f6
    @test convert(Q0f7, Rational{Int16}(-3//4)) === -0.75Q0f7
    @test convert(Q0f7, Rational{UInt8}(3//4)) === 0.75Q0f7
    @test_throws ArgumentError convert(Q0f7, typemax(Rational{Int8}))

    @test convert(Q0f7, Base.TwicePrecision(0.5)) === 0.5Q0f7
    @test_throws ArgumentError convert(Q7f8, Base.TwicePrecision(0x80, 0x01))
    tp = Base.TwicePrecision(0xFFFFFFFFp-32, 0xFFFFFFFEp-64)
    @test convert(Q0f63, tp) === reinterpret(Q0f63, typemax(Int64))
end

@testset "bool conversions" begin
    @test convert(Bool, 0.0Q1f6) === false
    @test convert(Bool, 1.0Q1f6) === true
    @test_throws InexactError convert(Bool, 0.5Q1f6)
    @test_throws InexactError convert(Bool, -1Q1f6)
    @test_broken convert(Bool, Fixed{Int8,8}(0.2)) # TODO: remove this
end

@testset "integer conversions" begin
    @test convert(Int, Q1f6(1)) === 1
    @test convert(Integer, Q1f6(1)) === Int8(1)
    @test convert(UInt, 1Q1f6) === UInt(1)
    @test_throws InexactError convert(Integer, 0.5Q1f6)
    @test_throws InexactError convert(Int8, 256Q9f6)
end

@testset "rational conversions" begin
    @test convert(Rational, -0.75Q1f6) === Rational{Int8}(-3//4)
    @test convert(Rational, -0.75Q0f7) === Rational{Int16}(-3//4)
    @test convert(Rational{Int}, -0.75Q0f7) === Rational{Int}(-3//4)

    @test rationalize(-0.75Q3f4) === Rational{Int}(-3//4)
    @test rationalize(Int16, 0.81Q3f4) === Rational{Int16}(13//16)
    @test rationalize(-0.81Q3f4, tol=0.02) === Rational{Int}(-13//16)
    @test rationalize(Int8, -0.81Q3f4, tol=0.07) === Rational{Int8}(-3//4)
end

@testset "BigFloat conversions" begin
    @test convert(BigFloat, -0.75Q0f7)::BigFloat == big"-0.75"

    @test big(Q7f0) === BigFloat # !== BigInt
    @test big(0.75Q3f4)::BigFloat == big"0.75"
end

@testset "float/floattype" begin
    @test float(0.75Q3f4) === 0.75f0
    @test float(0.75Q19f12) === 0.75
    @test float(0.75Q7f24) === 0.75
    @test float(0.75Q10f53)::BigFloat == big"0.75"

    @testset "floattype($F)" for F in target(Fixed, :i8, :i16, :i32, :i64; ex = :heavy)
        @test typemax(F) <= maxintfloat(floattype(F))
    end
end

@testset "conversions to float" begin
    for T in (Float16, Float32, Float64)
        @test isa(convert(T, Q0f7(0.3)), T)
    end

    for Tf in (Float16, Float32, Float64)
        @testset "$Tf(::$F)" for F in target(Fixed, :i8, :i16)
            T, f = rawtype(F), nbitsfrac(F)
            float_err = 0.0
            for i = typemin(T):typemax(T)
                f_expected = Tf(i * BigFloat(2)^-f)
                f_actual = Tf(reinterpret(F, i))
                float_err += abs(f_actual - f_expected)
            end
            @test float_err == 0.0
        end
        @testset "$Tf(::$F)" for F in target(Fixed, :i32, :i64, :i128)
            T, f = rawtype(F), nbitsfrac(F)
            error_count = 0
            for i in vcat(typemin(T):(typemin(T)+0xFF),
                          -T(0xFF):T(0xFF),
                          (typemax(T)-0xFF):typemax(T))
                f_expected = Tf(i * BigFloat(2)^-f)
                isinf(f_expected) && break # for Float16() and Float32()
                f_actual = Tf(reinterpret(F, i))
                f_actual == f_expected && continue
                error_count += 1
            end
            @test error_count == 0
        end
    end
end

@testset "fractional fixed-point numbers" begin
    # test all-fractional fixed-point numbers (issue #104)
    for F in (Q0f7, Q0f15, Q0f31, Q0f63)
        tmax = typemax(F)
        tol = (tmax + BigFloat(1.0)) / bitwidth(F)
        r = range(-1, stop=BigFloat(tmax)-tol, length=50)
        @test all(x -> abs(F(x) - x) <= tol, r)
    end
end

@testset "type modulus" begin
    T = Fixed{Int8,7}
    for i = -1.0:0.1:typemax(T)
        @test i % T === T(i)
    end
    @test ( 1.5 % T).i == round(Int,  1.5*128) % Int8
    @test (-0.3 % T).i == round(Int, -0.3*128) % Int8

    T = Fixed{Int16,9}
    for i = -64.0:0.1:typemax(T)
        @test i % T === T(i)
    end
    @test ( 65.2 % T).i == round(Int,  65.2*512) % Int16
    @test (-67.2 % T).i == round(Int, -67.2*512) % Int16
end

@testset "rounding" begin
    for sym in (:i8, :i16, :i32, :i64)
        T = symbol_to_inttype(Fixed, sym)
        rs = vcat([ oneunit(T) << b - oneunit(T) for b = 0:bitwidth(T)-1],
                  [ oneunit(T) << b              for b = 1:bitwidth(T)-2],
                  [ oneunit(T) << b + oneunit(T) for b = 2:bitwidth(T)-2],
                  [-oneunit(T) << b - oneunit(T) for b = 2:bitwidth(T)-2],
                  [-oneunit(T) << b              for b = 1:bitwidth(T)-1],
                  [-oneunit(T) << b + oneunit(T) for b = 1:bitwidth(T)-1])
        @testset "rounding $F" for F in target(Fixed, sym)
            xs = (reinterpret(F, r) for r in rs)
            @test all(x -> trunc(x) == trunc(float(x)), xs)
            @test all(x -> floor(float(x)) < typemin(F) || floor(x) == floor(float(x)), xs)
            @test all(x ->  ceil(float(x)) > typemax(F) ||  ceil(x) ==  ceil(float(x)), xs)
            @test all(x -> round(float(x)) > typemax(F) || round(x) == round(float(x)), xs)
            @test all(x -> trunc(Int64, x) === trunc(Int64, float(x)), xs)
            @test all(x -> floor(Int64, x) === floor(Int64, float(x)), xs)
            @test all(x ->  ceil(Int64, x) ===  ceil(Int64, float(x)), xs)
            @test all(x -> round(Int64, x) === round(Int64, float(x)), xs)
        end
    end
    @testset "rounding Fixed{Int16,$f} with overflow" for f = 1:16 # TODO: drop 16
        F = Fixed{Int16,f}
        @test_throws ArgumentError ceil(typemax(F))
        if f == 16
            @test_throws ArgumentError ceil(eps(F))
        elseif f == 15
            @test_throws ArgumentError ceil(eps(F))
            @test_throws ArgumentError round(typemax(F))
            @test_throws ArgumentError round(F(0.5) + eps(F))
        else
            @test_throws ArgumentError ceil(typemin(F) - oneunit(F) + eps(F))
            @test_throws ArgumentError round(typemax(F))
            @test_throws ArgumentError round(typemax(F) - F(0.5) + eps(F))
        end
    end
    @testset "rounding mode" begin
        @test round(-1.5Q1f6, RoundNearest) === -2Q1f6
        @test round(-1.5Q1f6, RoundToZero) === -1Q1f6
        @test round(-1.5Q1f6, RoundUp) === -1Q1f6
        @test round(-1.5Q1f6, RoundDown) === -2Q1f6
        @test round(Int, -1.5Q1f6, RoundNearest) === -2
        @test round(Int, -1.5Q1f6, RoundToZero) === -1
        @test round(Int, -1.5Q1f6, RoundUp) === -1
        @test round(Int, -1.5Q1f6, RoundDown) === -2
    end
    @test_throws InexactError trunc(UInt, typemin(Q0f7))
    @test_throws InexactError floor(UInt, -eps(Q0f7))
end

@testset "approx" begin
    @testset "approx $F" for F in target(Fixed, :i8, :i16; ex = :light)
        xs = typemin(F):eps(F):typemax(F)-eps(F)
        @test all(x -> x ≈ x + eps(F), xs)
        @test all(x -> x + eps(F) ≈ x, xs)
        @test !any(x -> x - eps(F) ≈ x + eps(F), xs)
    end
end

@testset "clamp" begin
    @test clamp(0.5Q0f7, -0.8Q0f7,  0.8Q0f7) === 0.5Q0f7
    @test clamp(0.5Q0f7, 0.75Q0f7,  0.8Q0f7) === 0.75Q0f7
    @test clamp(0.5Q0f7, -0.8Q0f7, 0.25Q0f7) === 0.25Q0f7
    @test clamp(0.5,      -0.8Q0f7,  0.8Q0f7) === 0.5
    @test clamp(0.5f0,    0.75Q0f7,  0.8Q0f7) === 0.75f0
    @test clamp(0.5Q0f15, -0.8Q0f7, 0.25Q0f7) === 0.25Q0f15
    @test clamp(0.5Q0f7, -Inf, Inf) === 0.5
    @test clamp(0.5,     Q0f7) === 0.5Q0f7
    @test clamp(-1.5f0,  Q0f7) === -1.0Q0f7
    @test clamp(1.5Q1f6, Q0f7) === 0.992Q0f7
end

@testset "sign-related functions" begin
    @test_throws Exception signed(Q0f7)
    @test_throws Exception signed(0.5Q0f7)
    @test_throws Exception unsigned(Q0f7)
    @test_throws Exception unsigned(0.5Q0f7)
    @test copysign(0.5Q0f7, 0x1) === 0.5Q0f7
    @test copysign(0.5Q0f7, -1) === -0.5Q0f7
    @test flipsign(0.5Q0f7, 0x1) === 0.5Q0f7
    @test flipsign(0.5Q0f7, -1) === -0.5Q0f7
    @test_throws ArgumentError sign(0Q0f7)
    @test sign(0Q1f6) === 0Q1f6
    @test sign(0.5Q1f6) === 1Q1f6
    @test sign(-0.5Q1f6) === -1Q1f6
    @test signbit(0.5Q0f7) === false
    @test signbit(-0.5Q0f7) === true
end

@testset "bitwise" begin
    @test bswap(Q0f7(0.5)) === Q0f7(0.5)
    @test bswap(Q0f15(0.5)) === reinterpret(Q0f15, signed(0x0040))
end

@testset "predicates" begin
    @test isfinite(1Q7f8)
    @test !isnan(1Q7f8)
    @test !isinf(1Q7f8)

    @testset "isinteger" begin
        @testset "isinteger(::$F)" for F in target(Fixed, :i8, :i16)
            xs = typemin(F):eps(F):typemax(F)
            @test all(x -> isinteger(x) == isinteger(float(x)), xs)
        end
        @testset "isinteger(::$F)" for F in target(Fixed, :i32, :i64, :i128)
            fzero, fmax, fmin = zero(F), typemax(F), typemin(F)
            if nbitsfrac(F) == 0
                @test isinteger(fzero) & isinteger(fmax) & isinteger(fmin)
            else
                @test isinteger(fzero) & !isinteger(fmax) & isinteger(fmin)
            end
        end
        @testset "isinteger(::Fixed{Int8,8})" begin # TODO: remove this testset
            @test !isinteger(Fixed{Int8,8}(-0.5))
            @test isinteger(Fixed{Int8,8}(0.0))
            @test !isinteger(Fixed{Int8,8}(127/256))
        end
    end
end

@testset "unit range" begin
    @test length(Q1f6(-1):Q1f6(0)) == 2
    @test length(Q1f6(0):Q1f6(-1)) == 0
    @test collect(Q1f6(-1):Q1f6(0)) == Q1f6[-1, 0]
    @test length(Q6f1(-64):Q6f1(63)) == 128
    QIntW = Fixed{Int,bitwidth(Int)-1}
    @test length(QIntW(-1):QIntW(0)) == 2
    QInt1 = Fixed{Int,1}
    @test length(typemin(QInt1):typemax(QInt1)-oneunit(QInt1)) == typemax(Int)
    @test_throws OverflowError length(typemin(QInt1):typemax(QInt1))
    @test length(-127Q7f0:127Q7f0) == 255
    @test length(Q1f62(0):Q1f62(-2)) == 0
end

@testset "step range" begin
    r = typemin(Q0f7):eps(Q0f7):typemax(Q0f7)
    counter = 0
    for x in r
        counter += 1
    end
    @test counter == 256
    @test length(r) == 256
    QInt1 = Fixed{Int,1}
    @test length(QInt1(0):eps(QInt1):typemax(QInt1)-eps(QInt1)) == typemax(Int)
    @test Base.unsafe_length(typemin(QInt1):eps(QInt1):typemax(QInt1)-eps(QInt1)) == -1
    @test_throws OverflowError length(QInt1(-1):eps(QInt1):typemax(QInt1)-eps(QInt1))
end

@testset "reductions" begin
    a = Q0f7[0.75, 0.5]
    acmp = Float64(a[1]) + Float64(a[2])
    @test sum(a) == acmp
    @test sum(a, dims=1) == [acmp]

    F6 = Fixed{Int8,6}
    a = F6[1.2, 1.4]
    acmp = Float64(a[1])*Float64(a[2])
    @test prod(a) == acmp
    @test prod(a, dims=1) == [acmp]
end

@testset "reductions, Statistics" begin
    a = Q1f6[0.75, 0.5]
    af = FixedPointNumbers.Treduce.(a)
    @test mean(a) === mean(af)
    @test std(a)  === std(af)
    @test var(a)  === var(af)
    m = mean(a)
    @test stdm(a, m) === stdm(af, m)
    @test varm(a, m) === varm(af, m)
end

@testset "rand" begin
    @testset "rand(::$F)" for F in target(Fixed; ex = :thin)
        @test isa(rand(F), F)
        a = rand(F, (3, 5))
        @test ndims(a) == 2 && eltype(a) === F
        @test size(a) == (3,5)
    end
    @test rand(MersenneTwister(1234), Q0f7) === -0.156Q0f7
end

@testset "Promotion within Fixed" begin
    @test @inferred(promote(Q0f7(0.25), Q0f7(0.75))) ===
        (Q0f7(0.25), Q0f7(0.75))
    @test @inferred(promote(Fixed{Int16,3}(0.25), Fixed{Int8,3}(0.875))) ===
        (Fixed{Int16,3}(0.25), Fixed{Int16,3}(0.875))
    @test @inferred(promote(Fixed{Int8,6}(0.125), Fixed{Int8,4}(0.75))) ===
        (Fixed{Int16,6}(0.125), Fixed{Int16,6}(0.75))

    @test Fixed{Int16,15}(-1)   == Fixed{Int8,7}(-1)
    @test Fixed{Int16,15}(0.25) == Fixed{Int8,7}(0.25)
    @test Fixed{Int16,7}(-1)   == Fixed{Int8,7}(-1)
    @test Fixed{Int16,7}(0.25) == Fixed{Int8,7}(0.25)
    @test Fixed{Int16,15}(-1)   == Fixed{Int8,5}(-1)
    @test Fixed{Int16,15}(5/32) == Fixed{Int8,5}(5/32)
    @test Fixed{Int16,3}(-1)   == Fixed{Int8,5}(-1)
    @test Fixed{Int16,3}(0.25) == Fixed{Int8,5}(0.25)

    @test @inferred(promote_type(Q0f7, Float64)) === Float64
    @test @inferred(promote_type(Float32, Q7f24)) === Float32 # Float64 on v0.9 (#207)

    @test @inferred(promote_type(Q0f7, Int8)) === Q0f7 # Float32 on v0.9 (#207)
    @test @inferred(promote_type(Int128, Q7f24)) === Q7f24 # Float64 on v0.9 (#207)

    @test @inferred(promote_type(Q0f15, Rational{UInt8})) === Rational{UInt8}

    @test @inferred(promote_type(Q0f7, Float32, Int)) === Float32
    @test @inferred(promote_type(Q0f7, Int, Float32)) === Float32
    @test @inferred(promote_type(Int, Q0f7, Float32)) === Float32
    @test @inferred(promote_type(Int, Float32, Q0f7)) === Float32
    @test @inferred(promote_type(Float32, Int, Q0f7)) === Float32
    @test @inferred(promote_type(Float32, Q0f7, Int)) === Float32

    @test @inferred(promote_type(Q0f7,Q1f6,Q2f5,Q3f4,Q4f3,Q5f2)) == Fixed{Int128,7}

    @test @inferred(promote_type(Q0f7, N0f32)) === FixedPoint # Float64 on v0.9 (#207)
end

@testset "show" begin
    iob = IOBuffer()
    q0f7 = reinterpret(Q0f7, signed(0xaa))
    show(iob, q0f7)
    str = String(take!(iob))
    @test str == "-0.672Q0f7"
    @test eval(Meta.parse(str)) === q0f7

    q15f16 = reinterpret(Q15f16, signed(0xaaaaaaaa))
    show(iob, q15f16)
    str = String(take!(iob))
    @test str == "-21845.33334Q15f16"
    @test eval(Meta.parse(str)) === q15f16

    show(IOContext(iob, :compact=>true), q15f16)
    @test String(take!(iob)) == "-21845.3"

    show(IOContext(iob, :compact=>true, :typeinfo=>Q15f16), q15f16)
    @test String(take!(iob)) == "-21845.3"

    show(IOContext(iob, :compact=>true, :typeinfo=>Fixed), q15f16)
    @test String(take!(iob)) == "-21845.3"

    show(IOContext(iob, :typeinfo=>Q15f16), q15f16)
    @test String(take!(iob)) == "-21845.33334Q15f16" # TODO: Consider removing suffix (issue #188)

    show(IOContext(iob, :typeinfo=>Normed), q15f16)
    @test String(take!(iob)) == "-21845.33334Q15f16"

    show(iob, Fixed{Int128,64}(-1.2345e6))
    @test String(take!(iob)) == "Fixed{Int128,$(SP)64}(-1.2345e6)"

    # TODO: remove this test
    show(iob, reinterpret(Fixed{Int8,8}, signed(0xaa)))
    @test String(take!(iob)) == "Fixed{Int8,$(SP)8}(-0.336)"
end

@testset "summary" begin
    a = Q0f7[0.2, 0.4]
    aa = Fixed[0.2Q0f7 0.4Q0f15]

    if VERSION >= v"1.6.0-DEV.356"
        @test_broken summary(a) == "2-element Vector{Q0f7}"
        @test_broken summary(view(a, 1:2)) == "2-element view(::Vector{Q0f7}, 1:2) with eltype Q0f7"
        @test_broken summary(aa) == "1×2 Matrix{Fixed}"
    else
        @test summary(a) == "2-element Array{Q0f7,1} with eltype Fixed{Int8,7}"
        @test summary(view(a, 1:2)) == "2-element view(::Array{Q0f7,1}, 1:2) with eltype Fixed{Int8,7}"
        @test_broken summary(aa) == "1×2 Array{Fixed,2}"
    end
end
