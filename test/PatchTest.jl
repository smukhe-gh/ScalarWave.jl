#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

# testing on a single patch
function check_calcPatch(N::Int, bnd1::Function, bnd2::Function)::Float64
    fpatch = Float64[bnd1(x) + bnd2(y) for x in chebgrid(N), y in chebgrid(N)]
    patch  = Patch([1,1], fpatch) 
    @test patch.value == fpatch
    @test patch.loc   == [1,1]
    @test getPB(patch, :R).value == fpatch[end, :]
    @test getPB(patch, :C).value == fpatch[:, end]
    bndx   = Boundary(:R, Float64[bnd1(x) for x in chebgrid(N)])
    bndy   = Boundary(:C, Float64[bnd2(y) for y in chebgrid(N)])
    npatch = calcPatch([1,1], bndx, bndy, operator(N,1)).value
    return L2norm(npatch, fpatch, chebweights(N))
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
function check_caclPatch_multipatch(N::Int, M::Int, loc::Array{Int,1}, bnd1::Function, bnd2::Function)::Float64
    #FIXME: flipping loc makes it work
    fpatch = Float64[bnd1(xg) + bnd2(yg) for xg in chebgrid(N, M, loc[2]), yg in chebgrid(N, M, loc[1])]
    bndx   = Boundary(:R, Float64[bnd1(xg) + bnd2(coordtransL2G(M, loc[2], 1.0)) for xg in chebgrid(N, M, loc[1])])
    bndy   = Boundary(:C, Float64[bnd1(coordtransL2G(M, loc[1], 1.0)) + bnd2(yg) for yg in chebgrid(N, M, loc[2])])
    npatch = calcPatch(loc, bndx, bndy, operator(N,1)).value
    return  L2norm(npatch, fpatch, (1/M)*chebweights(N))
end

@test check_caclPatch_multipatch(12, 2, [1,1], x->sin(pi*x), y->sin(pi*y)) < 1e-14
@test check_caclPatch_multipatch(12, 2, [2,2], x->sin(pi*x), y->sin(pi*y)) < 1e-14
@test check_caclPatch_multipatch(12, 2, [1,2], x->sin(pi*x), y->sin(pi*y)) < 1e-14
@test check_caclPatch_multipatch(12, 2, [2,1], x->sin(pi*x), y->sin(pi*y)) < 1e-14
