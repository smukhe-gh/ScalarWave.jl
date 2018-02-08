#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

@test 1==1

# testing p-convergence
for n in 2:20
    @show n, pconvergence(n, 2)
end

# testing h-convergence
for m in 2:20
    @show m, hconvergence(4,m)
end
