#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

module ScalarWave
using Luxor, Colors, FFTW, Roots, LinearAlgebra, PyPlot
using WriteVTK

export Grid, distribute
export Null, Spacelike, U, V,
       _uu, _dd, _u, _d, _udd

export Manifold, Space, ProductSpace, 
       Galerkin, Cardinal, 
       Chebyshev, GaussLobatto, Taylor, spacetype, Tag
export Field, Boundary, Operator, ProductSpaceOperator, IntegrationOperator,
       Metric, Derivative, Christoffel, Ricci,
       inversemetrictransform, metricdet
export order, dim, boundary, solve, ⦼, ⊙, shape, delta 
export collocation, derivative, integral,
       derivativetransform,
       basistransform, prefactor

export cheb, chebx, chebd, chebw, eye
export drawpatch, setcolormap, drawtensorfield
export find_t_of_UV, find_r_of_UV, find_U_of_tr, find_V_of_tr
export refine, driver, conductor 
export plot, pcolormesh, contourf, contour, levels

export F, J, Fvec, Jvec, Bvec, Svec, Sshape
export norm, Newton
export writevtk
export incomingboundary, outgoingboundary

include("datatypes/AbstractTypes.jl")
include("datatypes/Datatypes.jl")
include("dataTypes/MetricDataTypes.jl")

include("spectral/spaces/1DspaceNEW.jl")
include("spectral/spaces/2DspaceNEW.jl")

include("spectral/spaces/1Dspace.jl")
include("spectral/spaces/2Dspace.jl")
include("spectral/spaces/3Dspace.jl")
include("spectral/spaces/MetricSpaces.jl")
include("spectral/basis/TaylorGrid.jl")
include("spectral/basis/GaussLobattoGrid.jl")
include("spectral/basis/SpecBasis.jl")
include("spectral/basis/BasisTransformation.jl")
include("spectral/basis/Evaluate.jl")

include("amr/AMR.jl")
include("geometry/MetricFunctions.jl")
include("geometry/CoordinateTransform.jl")
include("spacetimes/schwarzschild/DoubleNullCoordinates.jl")
include("visualization/PyPlot.jl")
include("utilities/MathFunctions.jl")

include("einstein/FieldEquations.jl")
include("einstein/NonLinSolver.jl")
include("utilities/Utilities.jl")

end 
