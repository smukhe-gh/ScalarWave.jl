#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

module ScalarWave
using QuadGK

export Patch, Boundary
export cheb, chebx, chebd, chebw,
       pconvergence, hconvergence,
	   delta, coordtrans, reshapeA, 
       reshapeB, shapeA, shapeB, vandermonde, 
       prolongation1D, restriction1D,
       operator, getIC,
       getPB, pushforward, pullback, 
       distribute

include("SpecCalculus.jl")
include("Utilities.jl")
include("Physics.jl")
include("Patch.jl")
include("Grid.jl")

end 
