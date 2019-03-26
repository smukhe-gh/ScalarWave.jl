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

libraries = ["3DSpace"]
libraries = ["NLSolve"]
libraries = ["SpaceLikeHyperSurface"]
libraries = ["SkewedCoordinates"]
libraries = ["SummationByParts"]
libraries = ["SpecBasisCopy"]

for file in libraries
    @info "Testing $file"
    include("$(file)Test.jl")
end
