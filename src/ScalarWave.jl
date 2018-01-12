#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

module ScalarWave

export chebx, chebd, chebw,
	   operator, initializeB
#	   extractBC, setBC!,
#	   distribute,
#	   coordtrans, computeB

include("SpecCalculus.jl")
include("Physics.jl")
#include("Patch.jl")
#include("Grid.jl")
#include("Utilities.jl")

end 
