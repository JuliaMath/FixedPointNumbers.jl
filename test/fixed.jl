using Base.Test
using FixedPointNumbers
using Compat
import Compat.String

function test_op{F,T}(fun::F, ::Type{T}, fx, fy, fxf, fyf, tol)
    # Make sure that the result is representable
    (typemin(T) <= fun(fxf, fyf) <= typemax(T)) || return nothing
    @assert abs(fun(fx, fy) - convert(T, fun(fxf, fyf))) <= tol
    @assert abs(convert(Float64, fun(fx, fy)) - fun(fxf, fyf)) <= tol
end

function test_fixed{T}(::Type{T}, f)
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

@test isapprox(convert(Fixed{Int8,7}, 0.8), 0.797, atol=0.001)
@test isapprox(convert(Fixed{Int8,7}, 0.9), 0.898, atol=0.001)
@test_throws InexactError convert(Fixed{Int8, 7}, 0.999)
@test_throws InexactError convert(Fixed{Int8, 7}, 1.0)
@test_throws InexactError convert(Fixed{Int8, 7}, 1)
@test_throws InexactError convert(Fixed{Int8, 7}, 2)
@test_throws InexactError convert(Fixed{Int8, 7}, 128)

for (TI, f) in [(Int8, 8), (Int16, 8), (Int16, 10), (Int32, 16)]
    T = Fixed{TI,f}
    println("  Testing $T")
    test_fixed(T, f)
end

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

for T in [Fixed{Int8,7}, Fixed{Int16,8}, Fixed{Int16,10}]
    testapprox(T)  # defined in ufixed.jl
end

# reductions
F8 = Fixed{Int8,8}
a = F8[0.498, 0.1]
acmp = Float64(a[1]) + Float64(a[2])
@test sum(a) == acmp
@test sum(a, 1) == [acmp]

F6 = Fixed{Int8,6}
a = F6[1.2, 1.4]
acmp = Float64(a[1])*Float64(a[2])
@test prod(a) == acmp
@test prod(a, 1) == [acmp]

x = Fixed{Int8,8}(0.3)
for T in (Float16, Float32, Float64, BigFloat)
    y = convert(T, x)
    @test isa(y, T)
end

@test convert(Int, Q1f6(1)) === 1
@test convert(Integer, Q1f6(1)) === Int8(1)

# Floating-point conversions
@test isa(float(one(Fixed{Int8,6})),   Float32)
@test isa(float(one(Fixed{Int32,18})), Float64)
@test isa(float(one(Fixed{Int32,25})), Float64)

# Show
x = Fixed{Int32,5}(0.25)
iob = IOBuffer()
show(iob, x)
str = String(take!(iob))
@test str == "0.25Q26f5"
@test eval(parse(str)) == x

for T in (Fixed{Int8,8}, Fixed{Int16,8}, Fixed{Int16,10}, Fixed{Int32,16})
    a = rand(T)
    @test isa(a, T)
    a = rand(T, (3, 5))
    @test isa(a, Array{T,2})
    @test size(a) == (3,5)
end

# issue #79
@test realmin(Q11f4) == Q11f4(0.06)
