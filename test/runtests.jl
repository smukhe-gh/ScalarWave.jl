#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

using ScalarWave, Test

libraries = ["BasisTransformation",
             "MinkowskiDistorted",
             "Radial",
             "Ricci",
             "Schwarzschild"]

libraries = ["NonLinSolver"]
libraries = ["BasisTransformation"]
libraries = ["Evaluate"]

for file in libraries
    @info "Testing $file"
    include("$(file)Test.jl")
end
