#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

using ScalarWave
using Base.Test, PyPlot

libraries = ["MetricFunctions"]
libraries = ["Ricci"]
libraries = ["CoordinateTransform"]
libraries = ["Schwarzschild"]
libraries = ["SchwarzschildReggeWheeler"]
libraries = ["RicciNull"]
libraries = ["SchwarzschildEigenValues"]

for file in libraries
    info("Testing $file")
    include("$(file)Test.jl")
end
