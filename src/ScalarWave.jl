#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

module ScalarWave
using QuadGK

export Patch, Boundary
export cheb, chebx, chebd, chebw, chebgrid, chebweights, vandermonde,
	   delta, coordtransL2G, coordtransG2L, jacobian, 
       reshapeA, reshapeB, shapeA, shapeB,
       L2norm, L1norm, LInfnorm, savegrid, 
       getPB, calcPatch,  extractPatchCoeffs,
       interpolatePatch,
       operator, getIC,
       distribute

include("Types.jl")
include("SpecCalculus.jl")
include("Utilities.jl")
include("Patch.jl")
include("Grid.jl")
include("Physics.jl")
include("Dispatch.jl")
include("Visualization.jl")

end 
