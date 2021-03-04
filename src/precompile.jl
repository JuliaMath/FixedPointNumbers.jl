macro warnpcfail(ex::Expr)
    modl = __module__
    file = __source__.file === nothing ? "?" : String(__source__.file)
    line = __source__.line
    quote
        $(esc(ex)) || @warn """precompile directive
     $($(Expr(:quote, ex)))
 failed. Please report an issue in $($modl) (after checking for duplicates) or remove this directive.""" _file=$file _line=$line
    end
end

function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    normedtypes = (N0f8, N0f16)                      # precompiled Normed types
    realtypes = (Float16, Float32, Float64, Int)     # types for mixed Normed/Real operations
    for T in normedtypes
        for f in (+, -, abs, eps, rand)       # unary operations
            @warnpcfail precompile(Tuple{typeof(f),T})
        end
        @warnpcfail precompile(Tuple{typeof(rand),T,Tuple{Int}})
        @warnpcfail precompile(Tuple{typeof(rand),T,Tuple{Int,Int}})
        for f in (trunc, floor, ceil, round)  # rounding operations
            @warnpcfail precompile(Tuple{typeof(f),T})
            @warnpcfail precompile(Tuple{typeof(f),Type{Int},T})
        end
        for f in (+, -, *, /, <, <=, ==)      # binary operations
            @warnpcfail precompile(Tuple{typeof(f),T,T})
            for S in realtypes
                @warnpcfail precompile(Tuple{typeof(f),T,S})
                @warnpcfail precompile(Tuple{typeof(f),S,T})
            end
        end
        # conversions
        for S in realtypes
            @warnpcfail precompile(Tuple{Type{T},S})
            @warnpcfail precompile(Tuple{Type{S},T})
            @warnpcfail precompile(Tuple{typeof(convert),Type{T},S})
            @warnpcfail precompile(Tuple{typeof(convert),Type{S},T})
        end
    end
end
