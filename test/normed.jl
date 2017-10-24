using FixedPointNumbers
using Base.Test

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

UF2 = (Normed{UInt32,16}, Normed{UInt64,3}, Normed{UInt64,51}, Normed{UInt128,7}, Normed{UInt128,51})

for T in (FixedPointNumbers.UF..., UF2...)
    @test zero(T) == 0
    @test one(T) == 1
    @test one(T) * one(T) == one(T)
    @test typemin(T) == 0
    @test realmin(T) == eps(T)
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

x = N0f8(0.5)
@test convert(N0f8, x) === x
@test isfinite(x) == true
@test isnan(x) == false
@test isinf(x) == false

@test convert(N0f8,  1.1/typemax(UInt8)) == eps(N0f8)
@test convert(N6f10, 1.1/typemax(UInt16)*64) == eps(N6f10)
@test convert(N4f12, 1.1/typemax(UInt16)*16) == eps(N4f12)
@test convert(N2f14, 1.1/typemax(UInt16)*4)  == eps(N2f14)
@test convert(N0f16, 1.1/typemax(UInt16))    == eps(N0f16)
@test convert(Normed{UInt32,16}, 1.1/typemax(UInt32)*2^16) == eps(Normed{UInt32,16})
@test convert(Normed{UInt64,3},  1.1/typemax(UInt64)*UInt64(2)^61)  == eps(Normed{UInt64,3})
@test convert(Normed{UInt128,7}, 1.1/typemax(UInt128)*UInt128(2)^121) == eps(Normed{UInt128,7})

@test convert(N0f8,  1.1f0/typemax(UInt8)) == eps(N0f8)

@test convert(Float64, eps(N0f8)) == 1/typemax(UInt8)
@test convert(Float32, eps(N0f8)) == 1.0f0/typemax(UInt8)
@test convert(BigFloat, eps(N0f8)) == BigFloat(1)/typemax(UInt8)
for T in (FixedPointNumbers.UF..., UF2...)
    @test convert(Bool, zero(T)) == false
    @test convert(Bool, one(T))  == true
    @test convert(Bool, convert(T, 0.2)) == true
    @test convert(Int, one(T)) == 1
    @test convert(Integer, one(T)) == 1
    @test convert(Rational, one(T)) == 1
end
@test convert(Rational, convert(N0f8, 0.5)) == 0x80//0xff
@test convert(N0f16, one(N0f8)) === one(N0f16)
@test convert(N0f16, N0f8(0.5)).i === 0x8080
@test convert(Normed{UInt16,7}, Normed{UInt8,7}(0.504)) === Normed{UInt16,7}(0.504)

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

x = N0f8(0b01010001, 0)
@test ~x == N0f8(0b10101110, 0)
@test -x == reinterpret(N0f8, 0xaf)

@test isa(float(one(Normed{UInt8,7})),   Float32)
@test isa(float(one(Normed{UInt32,18})), Float64)
@test isa(float(one(Normed{UInt32,25})), Float64)

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

function testtrunc(inc::T) where {T}
    incf = convert(Float64, inc)
    tm = reinterpret(typemax(T))/reinterpret(one(T))
    x = zero(T)
    for i = 0 : min(1e6, reinterpret(typemax(T))-1)
        xf = incf*i
        try
            @test typeof(trunc(x)) == T
            @test trunc(x) == trunc(xf)
            @test typeof(round(x)) == T
            @test round(x) == round(xf)
            cxf = ceil(xf)
            if cxf < tm
                @test typeof(ceil(x)) == T
                @test ceil(x) == ceil(xf)
            end
            @test typeof(floor(x)) == T
            @test floor(x) == floor(xf)
            @test trunc(Int,x) == trunc(Int,xf)
            @test round(Int,x) == round(Int,xf)
            @test floor(Int,x) == floor(Int,xf)
            if cxf < tm
                @test ceil(Int,x) == ceil(Int,xf)
            end
        catch err
            println("Failed on x = ", x, ", xf = ", xf)
            rethrow(err)
        end
        x = convert(T, x+inc)
    end
end

for T in (FixedPointNumbers.UF..., UF2...)
    testtrunc(eps(T))
end

function testapprox(::Type{T}) where {T}
    for x = typemin(T):eps(T):typemax(T)-eps(T)
        y = x+eps(T)
        @test x ≈ y
        @test y ≈ x
        @test !(x ≈ y+eps(T))
    end
end
for T in FixedPointNumbers.UF
    testapprox(T)
end

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

r = reinterpret(N0f8, 0x01):reinterpret(N0f8, 0x01):reinterpret(N0f8, convert(UInt8, 48))
@test length(r) == 48

counter = 0
for x in N0f8(0):eps(N0f8):N0f8(1)
    local x
    counter += 1
end
@test counter == 256

# Promotion within Normed
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

# Show
x = reinterpret(N0f8, 0xaa)
iob = IOBuffer()
show(iob, x)
str = String(take!(iob))
@test str == "0.667N0f8"
@test eval(parse(str)) == x

# scaledual
function generic_scale!(C::AbstractArray, X::AbstractArray, s::Number)
    length(C) == length(X) || error("C must be the same length as X")
    for i = 1:length(X)
        @inbounds C[i] = X[i]*s
    end
    C
end

a = rand(UInt8, 10)
rfloat = similar(a, Float32)
rfixed = similar(rfloat)
af8 = reinterpret(N0f8, a)

b = 0.5
bd, eld = scaledual(b, af8[1])
@assert b*a[1] == bd*eld

b, ad = scaledual(0.5, a)
@test b == 0.5
@test ad == a
b, ad = scaledual(0.5, ad)
generic_scale!(rfloat, a, 0.5)
generic_scale!(rfixed, ad, b)
@test rfloat == rfixed

# reductions
a = N0f8[reinterpret(N0f8, 0xff), reinterpret(N0f8, 0xff)]
@test sum(a) == 2.0
@test sum(a, 1) == [2.0]

a = N2f14[3.2, 2.4]
acmp = Float64(a[1])*Float64(a[2])
@test prod(a) == acmp
@test prod(a, 1) == [acmp]

x = N0f8(0.3)
for T in (Float16, Float32, Float64, BigFloat)
    y = convert(T, x)
    @test isa(y, T)
end

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

# Overflow with Float16
@test N0f16(Float16(1.0)) === N0f16(1.0)
@test Float16(1.0) % N0f16 === N0f16(1.0)

if VERSION >= v"0.7.0-DEV.1790"
    a = N0f8[0.2, 0.4]
    @test summary(a) == "2-element Array{N0f8,1} with eltype FixedPointNumbers.Normed{UInt8,8}"
    @test summary(view(a, 1:2)) == "2-element view(::Array{N0f8,1}, 1:2) with eltype FixedPointNumbers.Normed{UInt8,8}"
end
