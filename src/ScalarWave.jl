#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

module ScalarWave
using Luxor, Colors, FFTW, Roots, LinearAlgebra, PyPlot

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
       basistransform, prefactor 

export cheb, chebx, chebd, chebw
export drawpatch, setcolormap, drawtensorfield
export find_t_of_UV, find_r_of_UV, find_U_of_tr, find_V_of_tr
export refine, coarsen
export plot, pcolormesh, contourf

include("utilities/datatypes/AbstractTypes.jl")
include("utilities/datatypes/Datatypes.jl")
include("utilities/dataTypes/MetricDataTypes.jl")

include("utilities/spectral/spaces/SingleSpaces.jl")
include("utilities/spectral/spaces/ProductSpaces.jl")
include("utilities/spectral/spaces/MetricSpaces.jl")
include("utilities/spectral/basis/TaylorGrid.jl")

include("utilities/spectral/basis/GaussLobattoGrid.jl")
include("utilities/spectral/basis/SpecBasis.jl")
include("utilities/spectral/basis/BasisTransformation.jl")
include("utilities/spectral/basis/Evaluate.jl")

include("utilities/AMR/AMR.jl")
include("geometry/MetricFunctions.jl")
include("geometry/CoordinateTransform.jl")
include("spacetimes/DoubleNullCoordinates.jl")
include("utilities/MathFunctions.jl")
include("utilities/visualization/PyPlot.jl")
end 
