#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Define the Abstract Types
#--------------------------------------------------------------------

abstract type Manifold{Tag} end
abstract type Space{Tag} <: Manifold{Tag} end
abstract type Cardinal{Tag} <: Space{Tag} end
abstract type Galerkin{Tag} <: Space{Tag} end

