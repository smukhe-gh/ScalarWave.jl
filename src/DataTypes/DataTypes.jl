#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Define concrete datatypes
#--------------------------------------------------------------------

struct U end
struct V end
struct UV end

struct uu end
struct dd end

struct Null end
struct Space end 

struct GaussLobatto{Tag ,N} <: Cardinal{Tag, N} end
struct Chebyshev{Tag ,N} <: Galerkin{Tag, N} end
struct Taylor{Tag ,N} <: Cardinal{Tag, N} end

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

struct Derivative{Tag, D}
    components::Array{T, 1} where {T}
end

struct Metric{Tag, D}
    components::Array{T, 1} where {T}
end

mutable struct Christoffel{Tag, D}
    components::Array{Field, 3}
end

mutable struct Ricci{Tag, D}
    components::Array{Field, 1}
end

# TODO: Refine patch
struct Patch
    loc::Array{Int,1}
    value::Array{Float64,2}
end
