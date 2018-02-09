#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

module ScalarWave
using QuadGK

export Patch, Boundary
export cheb, chebx, chebd, chebw,
       hconvergence, pconvergence,
	   delta, coordtrans, reshapeA, 
       reshapeB, shapeA, shapeB, vandermonde, chebgrid,
       operator, getIC,
       getPB, extractPatchCoeffs, interpolatePatch,
       distribute

include("SpecCalculus.jl")
include("Utilities.jl")
include("Physics.jl")
include("Patch.jl")
include("Grid.jl")

end 
