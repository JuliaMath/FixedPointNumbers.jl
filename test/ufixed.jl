using FixedPoint, Base.Test

# This is useful for testing
import FixedPoint: asraw
asraw(x) = x

@test asraw(0xa2uf8)  == 0xa2
@test asraw(0xa2uf10) == 0xa2
@test asraw(0xa2uf12) == 0xa2
@test asraw(0xa2uf14) == 0xa2
@test asraw(0xa2uf16) == 0xa2

for T in FixedPoint.UF
    @test zero(T) == 0
    @test one(T) == 1
    @test typemin(T) == 0
    @test realmin(T) == 0
    @test eps(zero(T)) == eps(typemax(T))
    @test sizeof(T) == 1 + (T != Ufixed8)
end
@test typemax(Ufixed8) == 1
@test typemax(Ufixed10) == typemax(Uint16)//(2^10-1)
@test typemax(Ufixed12) == typemax(Uint16)//(2^12-1)
@test typemax(Ufixed14) == typemax(Uint16)//(2^14-1)
@test typemax(Ufixed16) == 1
@test typemax(Ufixed10) == float64(typemax(Uint16))/(2^10-1)
@test typemax(Ufixed12) == float64(typemax(Uint16))/(2^12-1)
@test typemax(Ufixed14) == float64(typemax(Uint16))/(2^14-1)

@test convert(Ufixed8,  1.1/typemax(Uint8)) == eps(Ufixed8)
@test convert(Ufixed10, 1.1/typemax(Uint16)*64) == eps(Ufixed10)
@test convert(Ufixed12, 1.1/typemax(Uint16)*16) == eps(Ufixed12)
@test convert(Ufixed14, 1.1/typemax(Uint16)*4)  == eps(Ufixed14)
@test convert(Ufixed16, 1.1/typemax(Uint16))    == eps(Ufixed16)

@test convert(Ufixed8,  1.1f0/typemax(Uint8)) == eps(Ufixed8)

@test convert(Float64, eps(Ufixed8)) == 1/typemax(Uint8)
@test convert(Float32, eps(Ufixed8)) == 1.0f0/typemax(Uint8)
@test convert(BigFloat, eps(Ufixed8)) == BigFloat(1)/typemax(Uint8)
for T in FixedPoint.UF
    @test convert(Bool, zero(T)) == false
    @test convert(Bool, one(T))  == true
    @test convert(Bool, convert(T, 0.2)) == true
    @test convert(Int, one(T)) == 1
    @test convert(Rational, one(T)) == 1
end
@test convert(Rational, convert(Ufixed8, 0.5)) == 0x80//0xff

for T in FixedPoint.UF
    x = T(0x10,0)
    y = T(0x25,0)
    @test y > x
    @test y != x
    @test x+y == T(0x35,0)
    @test (x+y) - x == y
    @test (x-y) + y == x
end

function testtrunc{T}(inc::T)
    incf = float64(inc)
    tm = asraw(typemax(T))/asraw(one(T))
    x = zero(T)
    for i = 0:asraw(typemax(T))-1
        xf = incf*i
        try
            @test trunc(x) == trunc(xf)
            @test round(x) == round(xf)
            cxf = ceil(xf)
            if cxf < tm
                @test ceil(x) == ceil(xf)
            end
            @test floor(x) == floor(xf)
            @test itrunc(x) == itrunc(xf)
            @test iround(x) == iround(xf)
            @test ifloor(x) == ifloor(xf)
            if cxf < tm
                @test iceil(x) == iceil(xf)
            end
        catch err
            println("Failed on x = ", x, ", xf = ", xf)
            rethrow(err)
        end
        x += inc
    end
end

for T in FixedPoint.UF
    testtrunc(eps(T))
end
