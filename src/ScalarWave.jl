#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

module ScalarWave
using QuadGK, HCubature

export Patch, Boundary, LocalC, GlobalC
export cheb, chebx, chebd, chebw, chebgrid, chebweights, vandermonde, pseudovandermonde,
	   delta, coordtransL2G, coordtransG2L, jacobian, shapeH2L, shapeL2H, LInfnorm, L1norm, L2norm, 
       derivOP, boundaryOP,
       getPatchIC, getPatchBnd, calcPatch, extractPatchCoeffs, interpolatePatch, projectonPatchBnd, projectonPatch,
       distribute

include("Types.jl")
include("SpecCalculus.jl")
include("Utilities.jl")
include("Physics.jl")
include("Patch.jl")
include("Grid.jl")
include("Dispatch.jl")
include("Visualization.jl")

end 
