#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

using ScalarWave, Test

libraries = ["AxisOperator"]
libraries = ["AnalyticJacobian"]
libraries = ["MinkowskiEvolve"]
libraries = ["NonlinearSolver"]
libraries = ["Mix"]
libraries = ["MinkowskiEvolve"]

for file in libraries
    @info "Testing $file"
    include("$(file)Test.jl")
end
