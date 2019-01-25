#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

module ScalarWave
using Luxor, Colors, FFTW, Einsum, Roots, LinearAlgebra, Plots

export Grid, distribute
export Null, Spacelike, U, V,
       _uu, _dd, _u, _d, _udd

export Manifold, Space, ProductSpace, 
       Galerkin, Cardinal, 
       Chebyshev, GaussLobatto, Taylor, spacetype 
export Field, Boundary, Operator, ProductSpaceOperator, IntegrationOperator,
       Metric, Derivative, Christoffel, Ricci,
       inversemetrictransform, metricdet
export order, dim, boundary, solve, ⦼, shape, delta 
export collocation, derivative, integral,
       derivativetransform,
       basistransform, interpolate,  mapmetricindex, eye,
       L2Error, L2ErrorRelative
export cheb, chebx, chebd, chebw
export drawpatch, setcolormap, drawtensorfield
export find_t_of_UV, find_r_of_UV, find_U_of_tr, find_V_of_tr

include("DataTypes/AbstractTypes.jl")
include("DataTypes/DataTypes.jl")
include("DataTypes/MetricDataTypes.jl")
include("Spaces/SingleSpaces.jl")
include("Spaces/ProductSpaces.jl")
include("Spaces/MetricSpaces.jl")
include("Basis/TaylorGrid.jl")
include("Basis/GaussLobattoGrid.jl")
include("Basis/SpecBasis.jl")
include("Basis/BasisTransformation.jl")
include("Error.jl")
include("MetricFunctions.jl")
include("MathFunctions.jl")
include("CoordinateTransform.jl")
include("DoubleNullCoordinates.jl")
include("Visualization.jl")

end 
