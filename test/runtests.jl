#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

using ScalarWave
using Base.Test

modules = ["Grid",
           "SpecCalculus",
           "Utilities",
           "Physics",
           "Patch"]

for file in modules
    info("Testing $file")
    include("$(file)Test.jl")
end
