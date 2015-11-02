using FixedPointNumbers, Base.Test

@test (0xa2uf8)[] == 0xa2
@test (0xa2uf10)[] == 0xa2
@test (0xa2uf12)[] == 0xa2
@test (0xa2uf14)[] == 0xa2
@test (0xa2uf16)[] == 0xa2

@test reinterpret(UFixed8, 0xa2) == 0xa2uf8
@test reinterpret(UFixed10, 0x1fa2) == 0x1fa2uf10
@test reinterpret(UFixed12, 0x1fa2) == 0x1fa2uf12
@test reinterpret(UFixed14, 0x1fa2) == 0x1fa2uf14
@test reinterpret(UFixed16, 0x1fa2) == 0x1fa2uf16

@test ufixed8(1.0) == 0xffuf8
@test ufixed8(0.5) == 0x80uf8
@test ufixed14(1.0) == 0x3fffuf14
@test ufixed12([2]) == UFixed12[0x1ffeuf12]

for T in FixedPointNumbers.UF
    @test zero(T) == 0
    @test one(T) == 1
    @test one(T) * one(T) == one(T)
    @test typemin(T) == 0
    @test realmin(T) == 0
    @test eps(zero(T)) == eps(typemax(T))
    @test sizeof(T) == 1 + (T != UFixed8)
end
@test typemax(UFixed8) == 1
@test typemax(UFixed10) == typemax(UInt16)//(2^10-1)
@test typemax(UFixed12) == typemax(UInt16)//(2^12-1)
@test typemax(UFixed14) == typemax(UInt16)//(2^14-1)
@test typemax(UFixed16) == 1
@test typemax(UFixed10) == typemax(UInt16) // (2^10-1)
@test typemax(UFixed12) == typemax(UInt16) // (2^12-1)
@test typemax(UFixed14) == typemax(UInt16) // (2^14-1)

x = UFixed8(0.5)
@test isfinite(x) == true
@test isnan(x) == false
@test isinf(x) == false

@test convert(UFixed8,  1.1/typemax(UInt8)) == eps(UFixed8)
@test convert(UFixed10, 1.1/typemax(UInt16)*64) == eps(UFixed10)
@test convert(UFixed12, 1.1/typemax(UInt16)*16) == eps(UFixed12)
@test convert(UFixed14, 1.1/typemax(UInt16)*4)  == eps(UFixed14)
@test convert(UFixed16, 1.1/typemax(UInt16))    == eps(UFixed16)

@test convert(UFixed8,  1.1f0/typemax(UInt8)) == eps(UFixed8)

@test convert(Float64, eps(UFixed8)) == 1/typemax(UInt8)
@test convert(Float32, eps(UFixed8)) == 1.0f0/typemax(UInt8)
@test convert(BigFloat, eps(UFixed8)) == BigFloat(1)/typemax(UInt8)
for T in FixedPointNumbers.UF
    @test convert(Bool, zero(T)) == false
    @test convert(Bool, one(T))  == true
    @test convert(Bool, convert(T, 0.2)) == true
    @test convert(Int, one(T)) == 1
    @test convert(Rational, one(T)) == 1
end
@test convert(Rational, convert(UFixed8, 0.5)) == 0x80//0xff

x = UFixed8(0b01010001, 0)
@test ~x == UFixed8(0b10101110, 0)
@test -x == 0xafuf8

for T in FixedPointNumbers.UF
    x = T(0x10,0)
    y = T(0x25,0)
    fx = convert(Float32, x)
    fy = convert(Float32, y)
    @test y > x
    @test y != x
    @test typeof(x+y) == T
    @test typeof((x+y)-y) == T
    @test typeof(x*y) == T
    @test typeof(x/y) == T
    @test_approx_eq(x+y, T(0x35,0))
    @test_approx_eq((x+y)-x, fy)
    @test_approx_eq((x-y)+y, fx)
    @test_approx_eq(x*y, convert(T, fx*fy))
    @test_approx_eq(x/y, convert(T, fx/fy))
    @test_approx_eq(x^2, convert(T, fx^2))
    @test_approx_eq(x^2.1f0, fx^2.1f0)
    @test_approx_eq(x^2.1, convert(Float64, x)^2.1)
end

function testtrunc{T}(inc::T)
    incf = convert(Float64, inc)
    tm = typemax(T)[] / one(T)[]
    x = zero(T)
    for i = 0:typemax(T)[]-1
        xf = incf*i
        try
            @test trunc(x) == trunc(xf)
            @test round(x) == round(xf)
            cxf = ceil(xf)
            if cxf < tm
                @test ceil(x) == ceil(xf)
            end
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

for T in FixedPointNumbers.UF
    testtrunc(eps(T))
end

# Show
x = 0xaauf8
iob = IOBuffer()
show(iob, x)
str = takebuf_string(iob)
@test startswith(str, "UFixed8(")
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
