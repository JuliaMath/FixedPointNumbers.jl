using FixedPointNumbers, Compat.Test

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
        isinteger(x) && @show x
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
    for (TI, f) in [(Int8, 8), (Int16, 8), (Int16, 10), (Int32, 16)]
        T = Fixed{TI,f}
        println("  Testing $T")
        test_fixed(T, f)
    end
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
    for T in [Fixed{Int8,7}, Fixed{Int16,8}, Fixed{Int16,10}]
        testapprox(T)  # defined in ufixed.jl
    end
end

@testset "reductions" begin
    F8 = Fixed{Int8,8}
    a = F8[0.498, 0.1]
    acmp = convert(Float64, a[1]) + convert(Float64, a[2])
    @test sum(a) == acmp
    @test sum(a, dims=1) == [acmp]

    F6 = Fixed{Int8,6}
    a = F6[1.2, 1.4]
    acmp = convert(Float64, a[1])*convert(Float64, a[2])
    @test prod(a) == acmp
    @test prod(a, dims=1) == [acmp]
end

@testset "convert result type" begin
    x = Fixed{Int8,8}(0.3)
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

@testset "Show" begin
    x = Fixed{Int32,5}(0.25)
    iob = IOBuffer()
    show(iob, x)
    str = String(take!(iob))
    @test str == "0.25Q26f5"
    @test eval(Meta.parse(str)) == x
end

@testset "rand" begin
    for T in (Fixed{Int8,8}, Fixed{Int16,8}, Fixed{Int16,10}, Fixed{Int32,16})
        a = rand(T)
        @test isa(a, T)
        a = rand(T, (3, 5))
        @test ndims(a) == 2 && eltype(a) == T
        @test size(a) == (3,5)
    end
end

@testset "realmin" begin
    # issue #79
    @test realmin(Q11f4) == Q11f4(0.06)
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
        tol = (tmax + BigFloat(1.0)) / (sizeof(T) * 8)
        for x in range(-1, stop=convert(BigFloat, tmax)-tol, length=50)
            @test abs(Fixed{T, f}(x) - x) <= tol
        end
    end
end
