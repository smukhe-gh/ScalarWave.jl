#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

using ScalarWave, Test

libraries = ["Dispatch"]
libraries = ["BasisTransformation",
             "MinkowskiDistorted",
             "Schwarzschild"] 

libraries = ["Schwarzschild"]
libraries = ["Operator"]
libraries = ["Radial"]

for file in libraries
    @info "Testing $file"
    include("$(file)Test.jl")
end
