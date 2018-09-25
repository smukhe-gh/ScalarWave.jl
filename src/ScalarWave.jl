#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

module ScalarWave
using Luxor, Colors, FFTW

export Null, Spatial
export Manifold, Space, ProductSpace, 
       Galerkin, Cardinal, 
       Chebyshev, GaussLobatto, Taylor, 
       Field, Boundary, Operator, ProductSpaceOperator

export order, dim, len, identity, boundary, solve, ⦼, shape 
export collocation, derivative, 
       inversemetrictransform, inversemetricdet, derivativetransform,
       basistransform, interpolate

export Patch
export cheb, chebx, chebd, chebw, chebgrid,
	   delta, drawpatch

include("AbstractTypes.jl")
include("SingleSpaces.jl")
include("ProductSpaces.jl")
include("MathFunctions.jl")
include("CoordTransform.jl")
include("SpecCalculus.jl")
include("Visualization.jl")
include("Rationals.jl")
include("BasisTransformation.jl")

end 
