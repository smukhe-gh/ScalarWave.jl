#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

@test 1==1

@test pconvergence(20) < 1e-14
@test_broken hconvergence(4)  < 1e-14
