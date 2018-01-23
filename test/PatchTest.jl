#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

A = randn(5,5)
loc = [3,4]
patch = Patch(loc, A)
@test patch.loc   == loc
@test patch.value == A 

A[1,1] = 1.0
patch.value[1,1] == 1.0
@test patch.value == A

#bnd0 = randn(5)
#bnd1 = randn(5)
#patchboundary = Boundary(bnd0, bnd1)
#@test_broken patchboundary.row == bnd0
#@test_broken patchboundary.col == bnd1

#@test_broken getPB(patch, :R) == A[end, :]
#@test_broken getPB(patch, :C) == A[:, end]

#B = randn(5, 5)
#newpatch = Patch(loc, B)
#B[1, :] = bnd0
#B[:, 1] = bnd1
#setPB(newpatch, patchboundary)
#@test_broken newpatch.value == B

