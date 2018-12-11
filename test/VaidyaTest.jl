#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 12-2018
# Test functions for representing the Vaidya Metric
#--------------------------------------------------------------------

# basic tests of new functions to compute r(u,v) and f(u,v)
@test abs(r_of_m_exp(-40, 6, 1, 1) - 35.380170345966214) < 1e-14
@test abs(f_of_m_exp(-40, 6, 1, 1) - 0.9768937296256338) < 1e-14
rv = r_of_collapse(2.5, -3, 4, 1, 1)
@show Float64(abs(rv(2) - 3.019489599651289))
