#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Define the Abstract Types
#--------------------------------------------------------------------

abstract type Manifold{Tag} end
abstract type Space{Tag} <: Manifold{Tag} end
abstract type Galerkin{Tag, N} <: Space{Tag} end
abstract type Cardinal{Tag, N} <: Space{Tag} end

struct Patch
    loc::Array{Int,1}
    value::Array{Float64,2}
end
