using Base.Test
using FixedPointNumbers

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

            @test fx==fy || x!=y
            @test fx<fy  || x>=y
            @test fx<=fy || x>y

            for fun in [+, -, *, /]
                # Make sure that the result is representable
                if !(typemin(T) <= fun(fxf, fyf) <= typemax(T))
                    continue
                elseif (fun == /) && fy != 0
                    @test abs(fun(fx, fy) - convert(T, fun(fxf, fyf))) <= tol
                    @test abs(convert(Float64, fun(fx, fy)) - fun(fxf, fyf)) <= tol
                end
            end

            @test isequal(fx,fy) == isequal(hash(fx),hash(fy))
        end
    end
end

for (TI, f) in [(Int8, 8), (Int16, 8), (Int16, 10), (Int32, 16)]
    T = Fixed{TI,f}
    println("  Testing $T")
    test_fixed(T, f)
end