#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

module ScalarWave
using FastGaussQuadrature, Luxor, Colors, Roots

export Patch, Boundary, LocalC, GlobalC
export cheb, chebx, chebd, chebw, chebgrid, chebweights, vandermonde, pseudovandermonde,
	   delta, coordtransL2G, coordtransG2L, jacobian, shapeH2L, shapeL2H, 
       LInfnorm, L1norm, L2norm, array2dict, dict2array, savegrid, loadgrid, 
       derivOP, boundaryOP, RHS, modal2nodal, nodal2modal,
       getPatchIC, getPatchBnd, calcPatch, extractPatchCoeffs, interpolatePatch, 
       restrictmodes!, prolongatemodes, restrictOP, prolongateOP, restrictPatch, prolongatePatch,
       projectonPatchBndbyRestriction, projectonPatchbyRestriction, 
       distribute, sconv, showconv,
       drawmultipatch,
       fdistribute, fgetPatchBnd, fgetPatchIC, fRHS, fcalcPatch,
       PIdistribute, PIpotential,
       CompactifiedMinkowskidistribute

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
include("../beta/PIpotential.jl")
include("../beta/Compactification.jl")

end 
