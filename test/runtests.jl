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
             "Futures"]

smodule = ["Physics"]
smodule = ["Convergence"]
smodule = ["Visualization"]
smodule = ["PIpotential"]
smodule = ["Compactification"]

for file in smodule
    info("Testing $file")
    include("$(file)Test.jl")
end
