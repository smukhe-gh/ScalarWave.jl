#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

module ScalarWave
using QuadGK, ParallelAccelerator

export Patch, Boundary
export cheb, chebx, chebd, chebw,
       hconvergence, pconvergence,
	   delta, coordtrans, reshapeA, 
       reshapeB, shapeA, shapeB, vandermonde, chebgrid, 
       chebweights, operator, getIC,
       getPB, extractPatchCoeffs, interpolatePatch,
       distribute,
       plotgrid

include("Patch.jl")
include("SpecCalculus.jl")
include("Utilities.jl")
include("Physics.jl")
include("Grid.jl")
include("Visualization.jl")

end 
