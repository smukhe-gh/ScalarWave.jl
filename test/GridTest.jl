#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

@test 1==1

# testing h-convergence
for n in 2:20
    @show n, pconvergence(n, 2)
end
