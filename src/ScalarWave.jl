#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

module ScalarWave
using Luxor, Colors, JLD, PyPlot

export Manifold, Space, ProductSpace, 
       Galerkin, Cardinal, 
       Chebyshev, GaussLobatto, Taylor, 
       Field, Operator, ProductSpaceOperator

export order, dim, 
       derivative, identity, boundary, solve 

include("Spaces/AbstractTypes.jl")
include("Spaces/SingleSpaces.jl")
include("Spaces/ProductSpaces.jl")

export deriv, collocation
include("Rationals.jl")

export Patch, Boundary
export cheb, chebx, chebd, chebw, chebgrid, chebweights, vandermonde, pseudovandermonde,
	   delta, coordtransL2G, coordtransG2L, jacobian, shapeH2L, shapeL2H, 
       LInfnorm, L1norm, L2norm, array2dict, dict2array, savegrid, loadgrid, 
       derivOP, boundaryOP, RHS, modal2nodal, nodal2modal,
       getPatchIC, getPatchBnd, calcPatch, extractPatchCoeffs, interpolatePatch, 
       restrictmodes!, prolongatemodes, restrictOP, prolongateOP, restrictPatch, prolongatePatch,
       projectonPatchBndbyRestriction, projectonPatchbyRestriction, 
       distribute, drawmultipatch,
       fdistribute, fgetPatchBnd, fgetPatchIC, fRHS, fcalcPatch

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

end 
