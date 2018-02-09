#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

@test 1==1

#=
fr(x)    = sin(pi*x) 
fc(y)    = sin(pi*y)
sol(x,y) = sin(pi*x) + sin(pi*y)

errorvec = []
for n in 3:20
    L2errornorm = pconvergence(n, fr, fc, sol)
    @show n, L2errornorm
    append!(errorvec, L2errornorm)
end

# Hacky way to save the file
#using HDF5
#h5open("pconv.h5", "w") do file
#    g = g_create(file, "pconvergence") # create a group
#    g["n"]  = collect(3:20)              # create a scalar dataset inside the group
#    g["errorvec"] = errorvec  
#    attrs(g)["Description"] = "Test dataset" # an attribute
#end

temp = 1.0
for m in map(x -> 2^x, 2:9)
    error = hconvergence(1, m, fr, fc, sol)
    temp  = error
end
=#
