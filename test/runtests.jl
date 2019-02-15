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

libraries = ["Evaluate"]
libraries = ["BasisTransformation"]
libraries = ["AMR"]
libraries = ["PyPlot"]
libraries = ["NonLinSolver"]

for file in libraries
    @info "Testing $file"
    include("$(file)Test.jl")
end
