#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

using ScalarWave, Test

libraries = ["Collapse"]
libraries = ["MultipatchCollapse"]
libraries = ["Cowling"]

for file in libraries
    @info "Testing $file"
    include("$(file)Test.jl")
end
