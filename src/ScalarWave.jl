#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

module ScalarWave
using Luxor, Colors, FFTW, Einsum

export Null, Spatial, uu, dd, u, d
export Manifold, Space, ProductSpace, 
       Galerkin, Cardinal, 
       Chebyshev, GaussLobatto, Taylor, 
       Field, Boundary, Operator, ProductSpaceOperator,
       Metric, Derivative, Christoffel, Ricci,
       mapmetricindex

export order, dim, len, identity, boundary, solve, ⦼, shape 
export collocation, derivative, 
       inversemetrictransform, inversemetricdet, derivativetransform,
       basistransform, interpolate,
       metricinverse

export find_t_of_uv, find_r_of_uv

export Patch
export cheb, chebx, chebd, chebw, chebgrid,
	   delta, drawpatch

include("AbstractTypes.jl")
include("SingleSpaces.jl")
include("ProductSpaces.jl")
include("MetricFunctions.jl")
include("MathFunctions.jl")
include("CoordTransform.jl")
include("SpecCalculus.jl")
include("Visualization.jl")
include("Rationals.jl")
include("BasisTransformation.jl")
include("Coordinates.jl")

end 
