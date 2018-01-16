#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

A = randn(5,5)
@test_broken extractBC(A, 1) == A[end, :]
@test_broken extractBC(A, 0) == A[:, end]

B = zeros(5,5)
C = zeros(5,5)
B[1, :] = A[end, :]
B[:, 1] = A[:, end]
@test_broken setBC!(C, A[end, :], 1) == B[1,:]
@test_broken setBC!(C, A[:, end], 0) == B[:,1]
    
