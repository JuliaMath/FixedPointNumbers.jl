using Base.Test

function test_fixed{T}(::Type{T}, f)
    for x = -50:.01:50
        isinteger(x) && @show x
        fx = convert(T,x)
        @test convert(T,float64(fx)) == fx
        @test convert(T,float64(-fx)) == -fx
        @test float64(-fx) == -float64(fx)

        fxf = float64(fx)

        for y = -50:.01:50
            @show y
            fy = convert(T,y)
            fyf = float64(fy)

            @assert fx==fy || x!=y
            @assert fx<fy  || x>=y
            @assert fx<=fy || x>y

            @assert abs((fx+fy)-convert(T,fxf+fyf)) <= 2.0^-f
            @assert abs((fx-fy)-convert(T,fxf-fyf)) <= 2.0^-f
            @assert abs((fx*fy)-convert(T,fxf*fyf)) <= 2.0^-f
            if fy != 0
                @assert abs((fx/fy)-convert(T,fxf/fyf)) <= 2.0^-f
            end

            @assert abs(float64(fx+fy)-(fxf+fyf)) <= 2.0^-f
            @assert abs(float64(fx-fy)-(fxf-fyf)) <= 2.0^-f
            @assert abs(float64(fx*fy)-(fxf*fyf)) <= 2.0^-f
            if fy != 0
                @assert abs(float64(fx/fy)-(fxf/fyf)) <= 2.0^-f
            end
        end
    end
end

for f in [8, 10, 16]
    test_fixed(Fixed32{f}, f)
end
