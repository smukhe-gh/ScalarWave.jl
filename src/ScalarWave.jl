#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

module ScalarWave
using Luxor, Colors, PyPlot

include("Spaces/AbstractTypes.jl")
include("Spaces/SingleSpaces.jl")
include("Spaces/ProductSpaces.jl")
include("Rationals.jl")


export Manifold, Space, ProductSpace, 
       Galerkin, Cardinal, 
       Chebyshev, GaussLobatto, Taylor, 
       Field, Operator, ProductSpaceOperator

export order, dim, len, identity, boundary, solve 

export collocation, derivative

include("Types.jl")
include("SpecCalculus.jl")
include("Utilities.jl")
include("Physics.jl")
include("Operator.jl")
include("Patch.jl")
include("Grid.jl")
include("Projection.jl")
include("Dispatch.jl")
include("Visualization.jl")
include("Futures.jl")

export Patch, Boundary
export cheb, chebx, chebd, chebw, chebgrid, chebweights, vandermonde, pseudovandermonde,
	   delta, coordtransL2G, coordtransG2L, jacobian, shapeH2L, shapeL2H, 
       LInfnorm, L1norm, L2norm, array2dict, dict2array,
       derivOP, boundaryOP, RHS, modal2nodal, nodal2modal,
       getPatchIC, getPatchBnd, calcPatch, extractPatchCoeffs, interpolatePatch, 
       restrictmodes!, prolongatemodes, restrictOP, prolongateOP, restrictPatch, prolongatePatch,
       projectonPatchBndbyRestriction, projectonPatchbyRestriction, 
       distribute, drawmultipatch,
       fdistribute, fgetPatchBnd, fgetPatchIC, fRHS, fcalcPatch


end 
