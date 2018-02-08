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

func(x,y) = sin(pi*x) + sin(pi*y) 
for n in 2:20
    (N,M)     = (n,2) 
    loc       = [1,2]
    fx_global_exact  = Float64[func(i,j) for i in chebgrid(N), j in chebgrid(N)] 
    fx_global_exact_to_patch12 = prolongation2D(Patch([1,1], fx_global_exact), M, loc).value
   
    # what are other functions doing?
    fx_global_num = distribute(N, 1, x-> sin(pi*x), y->sin(pi*y))
    fx_global_num_to_patch12 = prolongation2D(fx_global_num[[1,1]], 2, loc).value 

    # check the prolongation for location [1,2]
    xf    = Float64[coordtrans(M, [chebx(i,N),chebx(1,N)], loc)[1] for i in 1:N+1]
    yf    = Float64[coordtrans(M, [chebx(1,N),chebx(j,N)], loc)[2] for j in 1:N+1]
    fx_patch12_exact = Float64[func(i,j) for i in xf, j in yf]
    #= 
    @show n, maximum(abs.(fx_global_exact  - fx_global_num[[1,1]].value))
    @show n, maximum(abs.(fx_patch12_exact - fx_global_exact_to_patch12))
    @show n, maximum(abs.(fx_patch12_exact - fx_global_num_to_patch12))
    @show n, maximum(abs.(fx_global_exact_to_patch12 - fx_global_num_to_patch12))
    println("--------------------------------------------------------------------------------")
    =#
end    

# checking restriction routines
fxpatches = distribute(2, 2, x-> sin(pi*x), y->sin(pi*y))
fxglobal  = restriction2D(fxpatches, 2)
@show size(fxglobal.value)
@show fxglobal.value









