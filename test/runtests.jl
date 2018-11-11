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

libraries = ["1D-Derivatives"]
libraries = ["2D-Derivatives "]
libraries = ["SchwarzschildMathematica"]
libraries = ["Integral"]

for file in libraries
    @info "Testing $file"
    include("$(file)Test.jl")
end
