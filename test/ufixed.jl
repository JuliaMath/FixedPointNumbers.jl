using FixedPointNumbers
using Base.Test
using Compat

@test reinterpret(0xa2uf8)  == 0xa2
@test reinterpret(0xa2uf10) == 0xa2
@test reinterpret(0xa2uf12) == 0xa2
@test reinterpret(0xa2uf14) == 0xa2
@test reinterpret(0xa2uf16) == 0xa2

@test reinterpret(UFixed8, 0xa2) == 0xa2uf8
@test reinterpret(UFixed10, 0x1fa2) == 0x1fa2uf10
@test reinterpret(UFixed12, 0x1fa2) == 0x1fa2uf12
@test reinterpret(UFixed14, 0x1fa2) == 0x1fa2uf14
@test reinterpret(UFixed16, 0x1fa2) == 0x1fa2uf16

@test UFixed8(1.0) == 0xffuf8
@test UFixed8(0.5) == 0x80uf8
@test UFixed14(1.0) == 0x3fffuf14
v = @compat UFixed12.([2])
@test v == UFixed12[0x1ffeuf12]
@test isa(v, Vector{UFixed12})

UF2 = (UFixed{UInt32,16}, UFixed{UInt64,3}, UFixed{UInt64,51}, UFixed{UInt128,7}, UFixed{UInt128,51})

for T in (FixedPointNumbers.UF..., UF2...)
    @test zero(T) == 0
    @test one(T) == 1
    @test one(T) * one(T) == one(T)
    @test typemin(T) == 0
    @test realmin(T) == 0
    @test eps(zero(T)) == eps(typemax(T))
    @test sizeof(T) == sizeof(FixedPointNumbers.rawtype(T))
end
@test typemax(UFixed8) == 1
@test typemax(UFixed10) == typemax(UInt16)//(2^10-1)
@test typemax(UFixed12) == typemax(UInt16)//(2^12-1)
@test typemax(UFixed14) == typemax(UInt16)//(2^14-1)
@test typemax(UFixed16) == 1
@test typemax(UFixed10) == typemax(UInt16) // (2^10-1)
@test typemax(UFixed12) == typemax(UInt16) // (2^12-1)
@test typemax(UFixed14) == typemax(UInt16) // (2^14-1)
@test typemax(UFixed{UInt32,16}) == typemax(UInt32) // (2^16-1)
@test typemax(UFixed{UInt64,3}) == typemax(UInt64) // (2^3-1)
@test typemax(UFixed{UInt128,7}) == typemax(UInt128) // (2^7-1)
@test typemax(UFixed{UInt128,100}) == typemax(UInt128) // (UInt128(2)^100-1)

# TODO: change back to InexactError when it allows message strings
@test_throws ArgumentError UFixed8(2)
@test_throws ArgumentError UFixed8(255)
@test_throws ArgumentError UFixed8(0xff)
@test_throws ArgumentError UFixed16(2)
@test_throws ArgumentError UFixed16(0xff)
@test_throws ArgumentError UFixed16(0xffff)
@test_throws ArgumentError convert(UFixed8,  typemax(UFixed10))
@test_throws ArgumentError convert(UFixed16, typemax(UFixed10))
@test_throws ArgumentError convert(UFixed{UInt128,100}, 10^9)
@test_throws ArgumentError convert(UFixed{UInt128,100}, 10.0^9)

x = UFixed8(0.5)
@test isfinite(x) == true
@test isnan(x) == false
@test isinf(x) == false

@test convert(UFixed8,  1.1/typemax(UInt8)) == eps(UFixed8)
@test convert(UFixed10, 1.1/typemax(UInt16)*64) == eps(UFixed10)
@test convert(UFixed12, 1.1/typemax(UInt16)*16) == eps(UFixed12)
@test convert(UFixed14, 1.1/typemax(UInt16)*4)  == eps(UFixed14)
@test convert(UFixed16, 1.1/typemax(UInt16))    == eps(UFixed16)
@test convert(UFixed{UInt32,16}, 1.1/typemax(UInt32)*2^16) == eps(UFixed{UInt32,16})
@test convert(UFixed{UInt64,3},  1.1/typemax(UInt64)*UInt64(2)^61)  == eps(UFixed{UInt64,3})
@test convert(UFixed{UInt128,7}, 1.1/typemax(UInt128)*UInt128(2)^121) == eps(UFixed{UInt128,7})

@test convert(UFixed8,  1.1f0/typemax(UInt8)) == eps(UFixed8)

@test convert(Float64, eps(UFixed8)) == 1/typemax(UInt8)
@test convert(Float32, eps(UFixed8)) == 1.0f0/typemax(UInt8)
@test convert(BigFloat, eps(UFixed8)) == BigFloat(1)/typemax(UInt8)
for T in (FixedPointNumbers.UF..., UF2...)
    @test convert(Bool, zero(T)) == false
    @test convert(Bool, one(T))  == true
    @test convert(Bool, convert(T, 0.2)) == true
    @test convert(Int, one(T)) == 1
    @test convert(Rational, one(T)) == 1
end
@test convert(Rational, convert(UFixed8, 0.5)) == 0x80//0xff
@test convert(UFixed16, one(UFixed8)) === one(UFixed16)
@test convert(UFixed16, UFixed8(0.5)).i === 0x8080
@test convert(UFixed{UInt16,7}, UFixed{UInt8,7}(0.504)) === UFixed{UInt16,7}(0.504)

@test  UFixed8(0.2) % UFixed8  === UFixed8(0.2)
@test UFixed14(1.2) % UFixed16 === UFixed16(0.20002)
@test UFixed14(1.2) % UFixed8  === UFixed8(0.196)

for i = 0.0:0.1:1.0
    @test i % UFixed8 === UFixed8(i)
end
@test ( 1.5 % UFixed8).i == round(Int,  1.5*255) % UInt8
@test (-0.3 % UFixed8).i == round(Int, -0.3*255) % UInt8

for i = 0.0:0.1:64.0
    @test i % UFixed10 === UFixed10(i)
end
@test (65.2 % UFixed10).i == round(Int, 65.2*1023) % UInt16
@test (-0.3 % UFixed10).i == round(Int, -0.3*1023) % UInt16

@test 1 % UFixed8 == 1
@test 2 % UFixed8 == UFixed8(0.996)

x = UFixed8(0b01010001, 0)
@test ~x == UFixed8(0b10101110, 0)
@test -x == 0xafuf8

@test isa(float(one(UFixed{UInt8,7})),   Float32)
@test isa(float(one(UFixed{UInt32,18})), Float64)
@test isa(float(one(UFixed{UInt32,25})), Float64)

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
    @test (x*y) ≈ convert(T, fx*fy)
    @test (x/y) ≈ convert(T, fx/fy)
    @test (x^2) ≈ convert(T, fx^2)
    @test (x^2.1f0) ≈ fx^2.1f0
    @test (x^2.1) ≈ convert(Float64, x)^2.1
end

function testtrunc{T}(inc::T)
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

function testapprox{T}(::Type{T})
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

@test !(UFixed8(0.5) < UFixed8(0.5))
@test UFixed8(0.5) <= UFixed8(0.5)

@test div(0x10uf8, 0x02uf8) == fld(0x10uf8, 0x02uf8) == 8
@test div(0x0fuf8, 0x02uf8) == fld(0x0fuf8, 0x02uf8) == 7
@test Base.fld1(0x10uf8, 0x02uf8) == 8
@test Base.fld1(0x0fuf8, 0x02uf8) == 8
@test mod(0x10uf8, 0x02uf8) == rem(0x10uf8, 0x02uf8) == 0
@test mod(0x0fuf8, 0x02uf8) == rem(0x0fuf8, 0x02uf8) == 0x01uf8
@test mod1(0x10uf8, 0x02uf8) == 0x02uf8
@test mod1(0x0fuf8, 0x02uf8) == 0x01uf8

r = 1uf8:1uf8:48uf8
@test length(r) == 48

counter = 0
for x in UFixed8(0):eps(UFixed8):UFixed8(1)
    counter += 1
end
@test counter == 256

# Promotion within UFixed
@test @inferred(promote(UFixed8(0.2), UFixed8(0.8))) ===
    (UFixed8(0.2), UFixed8(0.8))
@test @inferred(promote(UFixed{UInt16,3}(0.2), UFixed{UInt8,3}(0.86))) ===
    (UFixed{UInt16,3}(0.2), UFixed{UInt16,3}(0.86))
@test @inferred(promote(UFixed{UInt8,7}(0.197), UFixed{UInt8,4}(0.8))) ===
    (UFixed{UInt16,7}(0.197), UFixed{UInt16,7}(0.8))

@test UFixed{UInt16,16}(1)   == UFixed{UInt8,8}(1)
@test UFixed{UInt16,16}(0.2) == UFixed{UInt8,8}(0.2)
@test UFixed{UInt16,8}(1)   == UFixed{UInt8,8}(1)
@test UFixed{UInt16,8}(0.2) == UFixed{UInt8,8}(0.2)
@test UFixed{UInt16,16}(1)   == UFixed{UInt8,6}(1)
@test UFixed{UInt16,16}(0.20635) == UFixed{UInt8,6}(0.20635)
@test UFixed{UInt16,4}(1)   == UFixed{UInt8,6}(1)
@test UFixed{UInt16,4}(0.2) == UFixed{UInt8,6}(0.2)

@test promote_type(UFixed8,Float32,Int) == Float32
@test promote_type(UFixed8,Int,Float32) == Float32
@test promote_type(Int,UFixed8,Float32) == Float32
@test promote_type(Int,Float32,UFixed8) == Float32
@test promote_type(Float32,Int,UFixed8) == Float32
@test promote_type(Float32,UFixed8,Int) == Float32

# Show
x = 0xaauf8
iob = IOBuffer()
show(iob, x)
str = String(take!(iob))
@test str == "UFixed{UInt8,8}(0.667)"

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
af8 = reinterpret(UFixed8, a)

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
a = UFixed8[0xffuf8, 0xffuf8]
@test sum(a) == 2.0
@test sum(a, 1) == [2.0]

a = UFixed14[3.2, 2.4]
acmp = Float64(a[1])*Float64(a[2])
@test prod(a) == acmp
@test prod(a, 1) == [acmp]

x = UFixed8(0.3)
for T in (Float16, Float32, Float64, BigFloat)
    y = convert(T, x)
    @test isa(y, T)
end

for T in (UFixed{UInt8,8}, UFixed{UInt8,6},
          UFixed{UInt16,16}, UFixed{UInt16,14},
          UFixed{UInt32,32}, UFixed{UInt32,30},
          UFixed{UInt64,64}, UFixed{UInt64,62})
    a = rand(T)
    @test isa(a, T)
    a = rand(T, (3, 5))
    @test isa(a, Array{T,2})
    @test size(a) == (3,5)
end
