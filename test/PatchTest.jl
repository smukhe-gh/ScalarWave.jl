#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------
using QuadGK

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

function check_interpolation(N, M, loc)
    xg = chebgrid(2*N)
    xp = Float64[coordtrans(M, [x, 0], loc)[1] for x in xg]
    yp = Float64[coordtrans(M, [0, x], loc)[2] for x in xg]
    fp = Float64[sin(pi*x) + sin(pi*y) for x in xp, y in yp]

    fgrid  = distribute(N, M, x->sin(pi*x), y->sin(pi*y))
    fpatch = interpolatePatch(fgrid[loc[end:-1:1]], xg, xg) 
end

check_interpolation(14, 2, [1,1])
check_interpolation(14, 2, [2,1])
