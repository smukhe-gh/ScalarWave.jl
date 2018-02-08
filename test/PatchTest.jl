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

#-----------------------------------------------------------
# Testing restriction and prolongation for the 1D case
#-----------------------------------------------------------
# XXX: More debugging needed, since for fx = exp(-x^3/0.001)
# the error is not dropping with increasing number of modes.
N = 40
M = 1
xglobal   = Float64[chebx(i, N) for i in 1:N+1]
fxglobal  = Float64[sin(pi*x) for x in xglobal]
fxglobalR = restriction1D(prolongation1D(fxglobal, M), M)
@test maximum(abs.(fxglobal - fxglobalR)) < 1e-14

N = 4
M = 2
xglobal   = Float64[chebx(i, N) for i in 1:N+1]
fxglobal  = Float64[exp(-x^2/0.1) for x in xglobal]
fxglobalR = restriction1D(prolongation1D(fxglobal, M), M)
@test maximum(abs.(fxglobal - fxglobalR)) < 1e-14

N = 24
M = 3
xglobal   = Float64[chebx(i, N) for i in 1:N+1]
fxglobal  = Float64[exp(-x^2) for x in xglobal]
fxglobalR = restriction1D(prolongation1D(fxglobal, M), M)
@test maximum(abs.(fxglobal - fxglobalR)) < 1e-14

N = 60
M = 23
xglobal   = Float64[chebx(i, N) for i in 1:N+1]
fxglobal  = Float64[exp(-x^2/0.01) for x in xglobal]
fxglobalR = restriction1D(prolongation1D(fxglobal, M), M)
@test maximum(abs.(fxglobal - fxglobalR)) < 1e-14

#-----------------------------------------------------------
# test 2D prolongation and restriction routines
#-----------------------------------------------------------

M = 2
N = 32
xg = Float64[chebx(i, N) for i in 1:N+1]
fxglobal = Float64[sin(pi*i) + sin(pi*j) for i in xg, j in xg] 
patch = Patch([1,1], fxglobal)
fpatchdbase = prolongation2D(patch, M)
fxc12 = fpatchdbase[[1,2]].value

# check the prolongation for location [1,2]
xp12  = (xg - 1)/2
yp12  = (xg + 1)/2
fxp12 = Float64[sin(pi*i) + sin(pi*j) for i in xp12, j in yp12]
@test maximum(abs.(fxc12 - fxp12)) < 1e-14

# check the restriction routine
fxc  = restriction2D(fpatchdbase, M)
fxcv = fxc.value
#@show size(fxcv)
#@show size(fxglobal)
#@show maximum(abs.(fxcv - fxglobal))

