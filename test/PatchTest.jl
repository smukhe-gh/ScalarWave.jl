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

#--------------------------------------------------------------------
# Testing interpolation
#--------------------------------------------------------------------

# set chebyshev points
xcheb = Float64[chebx(i,3) for i in 1:4]

# now, construct xvec and yvec
xvec = (xcheb - 1)/2
yvec = (xcheb + 1)/2

# now, construct a function starting from coefficents with 3 modes
cfs  = rand(4,4)

# now, construct the vandermonde matrices for these
vndmx = vandermonde(3,xvec)
vndmy = vandermonde(3,yvec)

# now, construct the function 
func  = vndmx*cfs*vndmy'

# now, extract the coefficents and compare (works)
ecfs  = inv(vndmx)*func*inv(vndmy')
@test prod(ecfs - cfs .< 1e-13) == true

# now, recalculate the function values and compare
efunc = vndmx*ecfs*vndmy'
@test prod(efunc - func .< 1e-13) == true

# now, specify an arbitrary function and compute it's coefficents
sfunc  = Float64[x + y^3 + y^3*x^2 for x in xvec, y in yvec]
scfs   = inv(vndmx)*sfunc*inv(vndmy')
esfunc = vndmx*scfs*vndmy'
@show maximum(abs.(sfunc - esfunc))

# now, specify a Gaussian grid 
using QuadGK
xgauss, w = gauss(14)
xgauss = - xgauss 
xgvec  = (xgauss - 1)/2
ygvec  = (xgauss + 1)/2
vndmgx = vandermonde(3,xgvec)
vndmgy = vandermonde(3,ygvec)

# now, construct the function on the Gaussian grid
gsfunc = vndmgx*scfs*vndmgy'

# now, construct the exact function on Gaussian nodes
gefunc = Float64[x + y^3 + y^3*x^2 for x in xgvec, y in ygvec] 

# compare
@show maximum(abs.(gsfunc - gefunc))







