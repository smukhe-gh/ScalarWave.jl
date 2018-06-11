#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

module ScalarWave
using Luxor, Colors, JLD, SymEngine, Roots, PyPlot
export Patch, Boundary, D, Params, Grid, VarList

export cheb, chebx, chebd, chebw, chebgrid, chebweights, vandermonde, pseudovandermonde,
	   delta, coordtransL2G, coordtransG2L, jacobian, shapeH2L, shapeL2H, 
       LInfnorm, L1norm, L2norm, array2dict, dict2array, savegrid, loadgrid, 
       derivOP, boundaryOP, RHS, modal2nodal, nodal2modal,
       getPatchIC, getPatchBnd, calcPatch, extractPatchCoeffs, interpolatePatch, 
       restrictmodes!, prolongatemodes, restrictOP, prolongateOP, restrictPatch, prolongatePatch,
       projectonPatchBndbyRestriction, projectonPatchbyRestriction, 
       distribute, sconv, showconv, drawmultipatch,
       fdistribute, fgetPatchBnd, fgetPatchIC, fRHS, fcalcPatch,
       dict2struct, setgrid, setvarlist, 
       find_UV_of_TR, find_TR_of_UV, plotcoordgrid

include("Types.jl")
include("SpecCalculus.jl")
include("Utilities.jl")
include("Physics.jl")
include("Patch.jl")
include("Grid.jl")
include("Projection.jl")
include("Dispatch.jl")
include("Convergence.jl")
include("Visualization.jl")
include("Futures.jl")
include("Schwarzschild.jl")
include("Coordinates.jl")

end 
