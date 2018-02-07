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

# XXX: More debugging needed, since for fx = exp(-x^3/0.001)
#      the error is not dropping with increasing number of modes.
# Testing restriction and prolongation for the 1D case
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

#=
# resolution in the global patch
N = 64
M = 2

xglobal  = Float64[chebx(i, N) for i in 1:N+1]
fxglobal = Float64[x^2 for x in xglobal]
vndm     = vandermonde(N, xglobal)
xpatchl  = (xglobal+1)/2
xpatchr  = (xglobal-1)/2
pxl      = vandermonde(N,xpatchl)    
pxr      = vandermonde(N,xpatchr)

# define the prolongation and restriction operator
Pl = pxl*inv(vndm)
Pr = pxr*inv(vndm)
Rl = pinv(Pl)
Rr = pinv(Pr)

# from a single patch to two patches
fxl = Pl*fxglobal
fxr = Pr*fxglobal

#@show fxglobal
#@show fxl
#@show fxr

# from two patchs to a single patch
fxlglobal = Rl*fxl
fxrglobal = Rr*fxr

#@show fxlglobal
#@show fxrglobal

# Combine the vectors and project onto a single patch.
fxpatch = vcat(fxl, fxr)
P       = vcat(Pl, Pr)
@show size(fxpatch)
@show size(pinv(P))
fxptog = pinv(P)*fxpatch
@show fxptog
@show maximum(abs.(fxptog-fxglobal))
=#








