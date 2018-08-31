#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
# Rationals Test
#--------------------------------------------------------------------

N = 5
equispacedgrid = [collocation(Rational, k, N) for k in 1:N+1]
chebgausslobattogrid = [collocation(Float64, k, N) for k in 1:N+1]
@test chebd(2, 3, 10) == deriv(Float64, 2, 3, 10)
