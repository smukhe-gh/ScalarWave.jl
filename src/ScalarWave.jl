#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

module ScalarWave
using FastGaussQuadrature, Luxor, Colors, Roots
using SymEngine

export Patch, Boundary, LocalC, GlobalC
export cheb, chebx, chebd, chebw, chebgrid, chebweights, vandermonde, pseudovandermonde,
	   delta, coordtransL2G, coordtransG2L, jacobian, shapeH2L, shapeL2H, 
       LInfnorm, L1norm, L2norm, array2dict, dict2array, savegrid, loadgrid, 
       derivOP, boundaryOP, RHS, modal2nodal, nodal2modal,
       getPatchIC, getPatchBnd, calcPatch, extractPatchCoeffs, interpolatePatch, 
       restrictmodes!, prolongatemodes, restrictOP, prolongateOP, restrictPatch, prolongatePatch,
       projectonPatchBndbyRestriction, projectonPatchbyRestriction, 
       distribute, sconv, showconv, drawmultipatch,
       fdistribute, fgetPatchBnd, fgetPatchIC, fRHS, fcalcPatch
       
export Params, Grid, Metric 
export find_UV_from_TR, find_TR_from_UV, createmesh, creategrid, setmetric,
       computeaction, computexpansion

include("Types.jl")
include("SpecCalculus.jl")
include("Utilities.jl")
include("Physics.jl")
include("Patch.jl")
include("Grid.jl")
include("Projection.jl")
include("Dispatch.jl")
include("Convergence.jl")
include("Futures.jl")
include("Visualization.jl")

include("../beta/Auxillary.jl")
include("../beta/Schwarzschild.jl")

end 
