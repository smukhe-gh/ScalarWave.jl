#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

using ScalarWave, Test

libraries = ["InitialData"]
libraries = ["NonLinSolver"]
libraries = ["Performance"]

for file in libraries
    @info "Testing $file"
    include("$(file)Test.jl")
end
