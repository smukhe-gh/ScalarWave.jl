#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

@test 1==1

for n in 3:24
    @show n, pconvergence(n)
end

for m in 2:2:24
   @show m, hconvergence(12, m)
end
