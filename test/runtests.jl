#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

using ScalarWave
using Base.Test, PyPlot

libraries = ["MinkowskiDistortedRotation"]

for file in libraries
    tic()
    info("Testing $file")
    include("$(file)Test.jl")
    toc()
end
