#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

module ScalarWave
using QuadGK, JLD, ParallelAccelerator

export Patch, Boundary
export cheb, chebx, chebd, chebw, chebgrid, chebweights, vandermonde,
	   delta, coordtrans, jacobian, reshapeA, reshapeB, shapeA, shapeB,
       L2norm, L1norm, LInfnorm, savegrid, 
       getPB, calcPatch,  extractPatchCoeffs,
       operator, getIC,
       distribute,
       plotgrid

include("SpecCalculus.jl")
include("Utilities.jl")
include("Patch.jl")
include("Physics.jl")
include("Dispatch.jl")
include("Visualization.jl")

end 
