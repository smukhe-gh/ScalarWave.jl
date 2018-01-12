#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

using ScalarWave
using Base.Test

modules = ["SpecCalculus"]
println("Running tests:")
 
for file in modules
	println("- $(file)")
	include("$(file)Test.jl")
end