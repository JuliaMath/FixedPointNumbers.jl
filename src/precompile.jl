function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    normedtypes = (N0f8, N0f16)                      # precompiled Normed types
    realtypes = (Float16, Float32, Float64, Int)     # types for mixed Normed/Real operations
    for T in normedtypes
        for f in (+, -, abs, eps, rand)       # unary operations
            @assert precompile(Tuple{typeof(f),T})
        end
        @assert precompile(Tuple{typeof(rand),T,Tuple{Int}})
        @assert precompile(Tuple{typeof(rand),T,Tuple{Int,Int}})
        for f in (trunc, floor, ceil, round)  # rounding operations
            @assert precompile(Tuple{typeof(f),T})
            @assert precompile(Tuple{typeof(f),Type{Int},T})
        end
        for f in (+, -, *, /, <, <=, ==)      # binary operations
            @assert precompile(Tuple{typeof(f),T,T})
            for S in realtypes
                @assert precompile(Tuple{typeof(f),T,S})
                @assert precompile(Tuple{typeof(f),S,T})
            end
        end
        # conversions
        for S in realtypes
            @assert precompile(Tuple{Type{T},S})
            @assert precompile(Tuple{Type{S},T})
            @assert precompile(Tuple{typeof(convert),Type{T},S})
            @assert precompile(Tuple{typeof(convert),Type{S},T})
        end
    end
end
