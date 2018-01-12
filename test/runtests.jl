#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

using ScalarWave
using Base.Test

modules = ["SpecCalculus",
		   "Physics"]

for file in modules
	include("$(file)Test.jl")
end
