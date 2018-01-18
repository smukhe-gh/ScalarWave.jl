#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

A = randn(5,5)
@test extractBC(A, :R) == A[end, :]
@test extractBC(A, :C) == A[:, end]

M = zeros(5,5)
N = zeros(5,5)
setBC!(M, extractBC(A, :R), :R)
setBC!(N, extractBC(A, :C), :C)
@test M[1,:] == A[end, :]
@test N[:,1] == A[:, end]


