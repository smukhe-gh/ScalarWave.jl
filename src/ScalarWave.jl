#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

module ScalarWave
using Luxor, Colors

export Null, Spatial
export Manifold, Space, ProductSpace, 
       Galerkin, Cardinal, 
       Chebyshev, GaussLobatto, Taylor, 
       Field, Boundary, Operator, ProductSpaceOperator

export order, dim, len, identity, boundary, solve, ⦼, shape 
export collocation, derivative

export Patch
export cheb, chebx, chebd, chebw, chebgrid,
	   delta, drawpatch

include("Types.jl")
include("Spaces/AbstractTypes.jl")
include("Spaces/SingleSpaces.jl")
include("Spaces/ProductSpaces.jl")
include("Spaces/Mathfunctions.jl")
include("SpecCalculus.jl")
include("Visualization.jl")
include("Rationals.jl")

end 
