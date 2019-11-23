#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Define the Abstract Types
#--------------------------------------------------------------------

export Manifold, Space, ProductSpace, SingleSpaces,
       Field, Operator, U, V

abstract type Manifold{Tag} end
abstract type Space{Tag, N} <: Manifold{Tag} end

struct U end
struct V end

struct Field{S, D, T}
    space::S
    value::AbstractArray{T, D}
end

struct Operator{S, D, T}
    space::S
    value::AbstractArray{T, D}
end

struct ProductSpace{T1, T2}
    S1::T1
    S2::T2
end
