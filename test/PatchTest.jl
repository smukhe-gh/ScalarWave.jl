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

@test getPB(patch, :R).value == A[end, :]
@test getPB(patch, :C).value == A[:, end]

x = Array(linspace(1,-1,20))
@test size(interpolatePatch(patch, x, x).value) == (20, 20)

# test patch interpolation
xgrid = Float64[chebx(i, 10) for i in 1:11]
B = Float64[x^2 + y^3 for x in xgrid, y in xgrid]
npatch = Patch(loc, B)
@test isapprox(interpolatePatch(npatch, xgrid, xgrid).value, npatch.value, atol=15) 
