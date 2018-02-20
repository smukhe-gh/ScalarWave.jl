#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

# testing on a single patch
function check_calcPatch(N::Int, bnd1::Function, bnd2::Function)::Float64
    fpatch = Float64[bnd1(x) + bnd2(y) for x in chebgrid(N), y in chebgrid(N)]
    patch  = Patch([1,1], fpatch) 
    bndx   = Boundary(:R, Float64[bnd1(x) for x in chebgrid(N)])
    bndy   = Boundary(:C, Float64[bnd2(y) for y in chebgrid(N)])
    error  = calcPatch([1,1], bndx, bndy, operator(N,1)).value - fpatch
    @test patch.value == fpatch
    @test patch.loc   == [1,1]
    @test getPB(patch, :R).value == fpatch[end, :]
    @test getPB(patch, :C).value == fpatch[:, end]
    return  L2norm(calcPatch([1,1], bndx, bndy, operator(N,1)).value, fpatch, chebweights(N))
end

# these ones don't interpolate. You could have errors coming from Gaussian interpolation.
# they also don't project the boundaries. You could have errors due to aliasing.
@test check_calcPatch(12, x->sin(pi*x), y->sin(pi*y)) < 1e-14
@test check_calcPatch(12, x->exp(-x^2/0.001), y->exp(-y^2/0.001)) < 1e-14
@test check_calcPatch(12, x->x^2-1, y->y^2-1) < 1e-14

function check_extractPatchCoeffs(N::Int)::Float64
    coeffs  = randn(N+1,N+1)
    vndm    = vandermonde(N,chebgrid(N))
    fpatch  = vndm*coeffs*vndm'
    patch   = Patch([1,1], fpatch )
    ncoeffs = extractPatchCoeffs(patch)
    return L2norm(ncoeffs, coeffs, chebweights(N))
end
@test check_extractPatchCoeffs(12) < 1e-14
@test check_extractPatchCoeffs(24) < 1e-14

# testing on multi-patch systems
