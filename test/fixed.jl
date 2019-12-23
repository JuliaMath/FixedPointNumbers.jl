using FixedPointNumbers, Test
using FixedPointNumbers: bitwidth

function test_op(fun::F, ::Type{T}, fx, fy, fxf, fyf, tol) where {F,T}
    # Make sure that the result is representable
    (typemin(T) <= fun(fxf, fyf) <= typemax(T)) || return nothing
    @assert abs(fun(fx, fy) - convert(T, fun(fxf, fyf))) <= tol
    @assert abs(convert(Float64, fun(fx, fy)) - fun(fxf, fyf)) <= tol
end

function test_fixed(::Type{T}, f) where {T}
    values = [-10:0.01:10; -180:.01:-160; 160:.01:180]
    tol = 2.0^-f
    for x in values
        # Ignore values outside the representable range
        # typemin <, otherwise for -(-0.5)  > typemax
        if !(typemin(T) < x <= typemax(T))
            continue
        end
        fx = convert(T,x)
        @test convert(T,convert(Float64, fx)) == fx
        @test convert(T,convert(Float64, -fx)) == -fx
        @test convert(Float64, -fx) == -convert(Float64, fx)

        fxf = convert(Float64, fx)

        rx = convert(Rational{BigInt},fx)
        @assert isequal(fx,rx) == isequal(hash(fx),hash(rx))

        for y in values
            if !(typemin(T) < y <= typemax(T))
                continue
            end

            fy = convert(T,y)
            fyf = convert(Float64, fy)

            @assert fx==fy || x!=y
            @assert fx<fy  || x>=y
            @assert fx<=fy || x>y

            test_op(+, T, fx, fy, fxf, fyf, tol)
            test_op(-, T, fx, fy, fxf, fyf, tol)
            test_op(*, T, fx, fy, fxf, fyf, tol)
            fy != 0 && test_op(/, T, fx, fy, fxf, fyf, tol)

            @assert isequal(fx,fy) == isequal(hash(fx),hash(fy))
        end
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

@testset "reinterpret" begin
    @test reinterpret(Q0f7, signed(0xa2)) === -0.734375Q0f7
    @test reinterpret(Q5f10, signed(0x00a2)) === 0.158203125Q5f10

    @test reinterpret(reinterpret(Q0f7, signed(0xa2))) === signed(0xa2)
    @test reinterpret(reinterpret(Q5f10, signed(0x00a2))) === signed(0x00a2)

    @test reinterpret(Int8, 0.5Q0f7) === signed(0x40)
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

@testset "inexactness" begin
    @test_throws InexactError Q0f7(-2)
    # TODO: change back to InexactError when it allows message strings
    @test_throws ArgumentError one(Q0f15)
    @test_throws ArgumentError oneunit(Q0f31)
    @test_throws ArgumentError one(Fixed{Int8,8}) # TODO: remove this at end of its support
end

@testset "conversion" begin
    @test isapprox(convert(Fixed{Int8,7}, 0.8), 0.797, atol=0.001)
    @test isapprox(convert(Fixed{Int8,7}, 0.9), 0.898, atol=0.001)
    @test_throws InexactError convert(Fixed{Int8, 7}, 0.999)
    @test_throws InexactError convert(Fixed{Int8, 7}, 1.0)
    @test_throws InexactError convert(Fixed{Int8, 7}, 1)
    @test_throws InexactError convert(Fixed{Int8, 7}, 2)
    @test_throws InexactError convert(Fixed{Int8, 7}, 128)
end

@testset "test_fixed" begin
    for (TI, f) in [(Int8, 7), (Int16, 8), (Int16, 10), (Int32, 16)]
        T = Fixed{TI,f}
        # println("  Testing $T")
        test_fixed(T, f)
    end
end

@testset "rounding" begin
    for T in (Int8, Int16, Int32, Int64)
        rs = vcat([ oneunit(T) << b - oneunit(T) for b = 0:bitwidth(T)-1],
                  [ oneunit(T) << b              for b = 1:bitwidth(T)-2],
                  [ oneunit(T) << b + oneunit(T) for b = 2:bitwidth(T)-2],
                  [-oneunit(T) << b - oneunit(T) for b = 2:bitwidth(T)-2],
                  [-oneunit(T) << b              for b = 1:bitwidth(T)-1],
                  [-oneunit(T) << b + oneunit(T) for b = 1:bitwidth(T)-1])
        @testset "rounding Fixed{$T,$f}" for f = 0:bitwidth(T)-1
            F = Fixed{T,f}
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

@testset "modulus" begin
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

@testset "testapprox" begin
    @testset "approx $T" for T in [Fixed{Int8,7}, Fixed{Int16,8}, Fixed{Int16,10}]
        xs = typemin(T):eps(T):typemax(T)-eps(T)
        @test all(x -> x ≈ x + eps(T), xs)
        @test all(x -> x + eps(T) ≈ x, xs)
        @test !any(x -> x - eps(T) ≈ x + eps(T), xs)
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

@testset "convert result type" begin
    x = Fixed{Int8,7}(0.75)
    for T in (Float16, Float32, Float64, BigFloat)
        y = convert(T, x)
        @test isa(y, T)
    end
end

@testset "Integer conversions" begin
    @test convert(Int, Q1f6(1)) === 1
    @test convert(Integer, Q1f6(1)) === Int8(1)
end

@testset "Floating-point conversions" begin
    @test isa(float(one(Fixed{Int8,6})),   Float32)
    @test isa(float(one(Fixed{Int32,18})), Float64)
    @test isa(float(one(Fixed{Int32,25})), Float64)
end

@testset "predicates" begin
    @test isfinite(1Q7f8)
    @test !isnan(1Q7f8)
    @test !isinf(1Q7f8)

    @testset "isinteger" begin
        for T in (Int8, Int16)
            @testset "isinteger(::Fixed{$T,$f})" for f = 0:bitwidth(T)-1
                F = Fixed{T,f}
                xs = typemin(F):eps(F):typemax(F)
                @test all(x -> isinteger(x) == isinteger(float(x)), xs)
            end
        end
        for T in (Int32, Int64)
            @testset "isinteger(::Fixed{$T,$f})" for f = 0:bitwidth(T)-1
                F = Fixed{T,f}
                fzero, fmax, fmin = zero(F), typemax(F), typemin(F)
                if f == 0
                    @test isinteger(fzero) & isinteger(fmax) & isinteger(fmin)
                else
                    @test isinteger(fzero) & !isinteger(fmax) & isinteger(fmin)
                end
            end
        end
        @testset "isinteger(::Fixed{Int8,8})" begin # TODO: remove this testset
            @test !isinteger(Fixed{Int8,8}(-0.5))
            @test isinteger(Fixed{Int8,8}(0.0))
            @test !isinteger(Fixed{Int8,8}(127/256))
        end
    end
end

@testset "Show" begin
    x = Fixed{Int32,5}(0.25)
    iob = IOBuffer()
    show(iob, x)
    str = String(take!(iob))
    @test str == "0.25Q26f5"
    @test eval(Meta.parse(str)) == x
end

@testset "rand" begin
    for F in (Fixed{Int8,7}, Fixed{Int16,8}, Fixed{Int16,10}, Fixed{Int32,16})
        @test isa(rand(F), F)
        a = rand(F, (3, 5))
        @test ndims(a) == 2 && eltype(a) == F
        @test size(a) == (3,5)
    end
end

@testset "floatmin" begin
    # issue #79
    @test floatmin(Q11f4) == Q11f4(0.06)
end

@testset "Disambiguation constructors" begin
    @test_throws ArgumentError Fixed{Int32,16}('a')
    @test_throws InexactError  Fixed{Int32,16}(complex(1.0, 1.0))
    @test Fixed{Int32,16}(complex(1.0, 0.0))        == 1
    @test Fixed{Int32,16}(Base.TwicePrecision(1.0, 0.0)) == 1
end

@testset "fractional fixed-point numbers" begin
    # test all-fractional fixed-point numbers (issue #104)
    for (T, f) in ((Int8, 7),
                 (Int16, 15),
                 (Int32, 31),
                 (Int64, 63))
        tmax = typemax(Fixed{T, f})
        @test tmax == BigInt(typemax(T)) / BigInt(2)^f
        tol = (tmax + BigFloat(1.0)) / bitwidth(T)
        for x in range(-1, stop=BigFloat(tmax)-tol, length=50)
            @test abs(Fixed{T, f}(x) - x) <= tol
        end
    end
end

@testset "low-level arithmetic" begin
    @test bswap(Q0f7(0.5)) === Q0f7(0.5)
    @test bswap(Q0f15(0.5)) === reinterpret(Q0f15, signed(0x0040))
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

    @test promote_type(Q0f7,Float32,Int) == Float32
    @test promote_type(Q0f7,Int,Float32) == Float32
    @test promote_type(Int,Q0f7,Float32) == Float32
    @test promote_type(Int,Float32,Q0f7) == Float32
    @test promote_type(Float32,Int,Q0f7) == Float32
    @test promote_type(Float32,Q0f7,Int) == Float32
    @test promote_type(Q0f7,Q1f6,Q2f5,Q3f4,Q4f3,Q5f2) == Fixed{Int128,7}
end
