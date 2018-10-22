#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Define concrete datatypes
#--------------------------------------------------------------------

struct U  end
struct V  end
struct UV end

struct uu end
struct dd end
struct u  end
struct d  end
struct udd end

struct Null end
struct Spacelike end 

struct GaussLobatto{Tag ,N} <: Space{Tag} end
struct Chebyshev{Tag ,N} <: Space{Tag} end
struct Taylor{Tag ,N} <:Space{Tag} end

Cardinal{Tag, N} = Union{GaussLobatto{Tag, N}, Taylor{Tag, N}} 

struct Field{S, D, T}
    space::Type{S}
    value::Array{T, D}
end

struct Boundary{S, D, T}
    space::Type{S}
    value::Array{T, D}
end

struct Operator{S, D, T}
    space::Type{S}
    value::Array{T, D}
end

struct ProductSpace{S1<:Space, S2<:Space} end

struct ProductSpaceOperator{S, D, T}
    space::Type{S}
    value::Array{T, D}
end

# TODO: Refine patch to work with coordinate extents
struct Patch
    loc::Array{Int,1}
    value::Array{Float64,2}
end
