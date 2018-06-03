#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

using ScalarWave
using Base.Test

libraries = ["SpecCalculus",
             "Utilities",
             "Physics",
             "Patch",
             "Grid",
             "Projection",
             "Visualization",
             "Dispatch",
             "Futures",
             "Derivatives"]

smodule = ["Convergence"]
smodule = ["Schwarzschild"]

for file in libraries
    info("Testing $file")
    include("$(file)Test.jl")
end
