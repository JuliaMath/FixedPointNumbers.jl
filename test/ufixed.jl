using FixedPointNumbers, Base.Test

@test reinterpret(0xa2uf8)  == 0xa2
@test reinterpret(0xa2uf10) == 0xa2
@test reinterpret(0xa2uf12) == 0xa2
@test reinterpret(0xa2uf14) == 0xa2
@test reinterpret(0xa2uf16) == 0xa2

@test reinterpret(Ufixed8, 0xa2) == 0xa2uf8
@test reinterpret(Ufixed10, 0x1fa2) == 0x1fa2uf10
@test reinterpret(Ufixed12, 0x1fa2) == 0x1fa2uf12
@test reinterpret(Ufixed14, 0x1fa2) == 0x1fa2uf14
@test reinterpret(Ufixed16, 0x1fa2) == 0x1fa2uf16

for T in FixedPointNumbers.UF
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
for T in FixedPointNumbers.UF
    @test convert(Bool, zero(T)) == false
    @test convert(Bool, one(T))  == true
    @test convert(Bool, convert(T, 0.2)) == true
    @test convert(Int, one(T)) == 1
    @test convert(Rational, one(T)) == 1
end
@test convert(Rational, convert(Ufixed8, 0.5)) == 0x80//0xff

for T in FixedPointNumbers.UF
    x = T(0x10,0)
    y = T(0x25,0)
    @test y > x
    @test y != x
    @test x+y == T(0x35,0)
    @test (x+y) - x == y
    @test (x-y) + y == x
    @test x*y == float32(x)*float32(y)
    @test x/y == float32(x)/float32(y)
    @test x^2 == float32(x)^2
    @test x^2.1f0 == float32(x)^2.1f0
    @test x^2.1 == float64(x)^2.1
end

function testtrunc{T}(inc::T)
    incf = float64(inc)
    tm = reinterpret(typemax(T))/reinterpret(one(T))
    x = zero(T)
    for i = 0:reinterpret(typemax(T))-1
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

for T in FixedPointNumbers.UF
    testtrunc(eps(T))
end

# scaledual
function generic_scale!(C::AbstractArray, X::AbstractArray, s::Number)
    length(C) == length(X) || error("C must be the same length as X")
    for i = 1:length(X)
        @inbounds C[i] = X[i]*s
    end
    C
end

a = rand(Uint8, 10)
rfloat = similar(a, Float32)
rfixed = similar(rfloat)
af8 = reinterpret(Ufixed8, a)

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
