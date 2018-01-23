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

bnd0 = randn(5)
bnd1 = randn(5)
patchboundary = Boundary(bnd0, bnd1)
@test patchboundary.row == bnd0
@test patchboundary.col == bnd1

@test getPB(patch, :R) == A[end, :]
@test getPB(patch, :C) == A[:, end]

B = randn(5, 5)
newpatch = Patch(loc, B)
B[1, :] = bnd0
B[:, 1] = bnd1
setPB(newpatch, patchboundary)
@test newpatch.value == B

function modpatch(patch::Patch)::Patch
    patch.value = patch.value*3.0
    return patch
end

X = randn(6,6)
testpatch = Patch([0,0], X)
modpatch(testpatch)
@test testpatch.value == X*3.0
