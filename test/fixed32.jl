using Base.Test
using FixedPointNumbers

function test_fixed{T}(::Type{T}, f)
    values = [-10:0.01:10, -180:.01:-160, 160:.01:180]
    tol = 2.0^-f

    for x in values
        #isinteger(x) && @show x
        fx = convert(T,x)
        @test convert(T,float64(fx)) == fx
        @test convert(T,float64(-fx)) == -fx
        @test float64(-fx) == -float64(fx)

        fxf = float64(fx)

        for y in values
            fy = convert(T,y)
            fyf = float64(fy)

            @assert fx==fy || x!=y
            @assert fx<fy  || x>=y
            @assert fx<=fy || x>y

            @assert abs((fx+fy)-convert(T,fxf+fyf)) <= tol
            @assert abs((fx-fy)-convert(T,fxf-fyf)) <= tol
            @assert abs((fx*fy)-convert(T,fxf*fyf)) <= tol
            if fy != 0
                @assert abs((fx/fy)-convert(T,fxf/fyf)) <= tol
            end

            @assert abs(float64(fx+fy)-(fxf+fyf)) <= tol
            @assert abs(float64(fx-fy)-(fxf-fyf)) <= tol
            @assert abs(float64(fx*fy)-(fxf*fyf)) <= tol
            if fy != 0
                @assert abs(float64(fx/fy)-(fxf/fyf)) <= tol
            end
        end
    end
end

for f in [8, 10, 16]
    test_fixed(Fixed32{f}, f)
end
