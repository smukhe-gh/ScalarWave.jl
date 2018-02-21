#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

using ScalarWave
using Base.Test, QuadGK

modules = ["SpecCalculus",
           "Physics"]

for file in modules
    info("Testing $file")
    include("$(file)Test.jl")
end
