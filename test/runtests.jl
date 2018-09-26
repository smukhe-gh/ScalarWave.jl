#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

using ScalarWave
using Base.Test, PyPlot

libraries = ["Spaces"]
libraries = ["Wave"]
libraries = ["BasisTransformation"]
libraries = ["Schwarzschild"]
libraries = ["MetricFunctions"]

for file in libraries
    info("Testing $file")
    include("$(file)Test.jl")
end
