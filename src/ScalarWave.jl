#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

module ScalarWave

export chebx, chebd, chebw,
	   delta, coordtrans, reshapeA, reshapeB, shapeA, shapeB, computeRHS,
       operator, initializeRHS,
       extractBC, setBC!

include("SpecCalculus.jl")
include("Utilities.jl")
include("Physics.jl")
include("Patch.jl")

end 
