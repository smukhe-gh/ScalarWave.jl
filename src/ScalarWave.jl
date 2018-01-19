#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

module ScalarWave

export chebx, chebd, chebw,
	   delta, coordtrans, reshapeA, reshapeB, shapeA, shapeB, computeRHS,
       operator, operatorNBC, operatorNBCW, initializeRHS,
       extractBC, setBC!,
       fextractBC, fsetBC

include("SpecCalculus.jl")
include("Utilities.jl")
include("Physics.jl")
include("Patch.jl")
include("Futures.jl")
end 
