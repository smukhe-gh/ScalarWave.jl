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
             "Dispatch",
             "Futures",
             "Visualization"]

libraries = ["Coordinates"]
libraries = ["Schwarzschild"]

for file in libraries
    info("Testing $file")
    include("$(file)Test.jl")
end
