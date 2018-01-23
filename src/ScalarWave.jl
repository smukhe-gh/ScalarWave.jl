#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

module ScalarWave

export Patch, Boundary

export chebx, chebd, chebw,
	   delta, coordtrans, reshapeA, reshapeB, shapeA, shapeB, 
       operator, getIC,
       getPB, setPB

include("SpecCalculus.jl")
include("Utilities.jl")
include("Physics.jl")
include("Patch.jl")

end 
