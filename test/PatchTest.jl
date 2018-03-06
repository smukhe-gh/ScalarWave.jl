#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function testextractPatchCoeffs(Nx::Int, Ny::Int)::Float64
    coeffs  = randn(Nx+1,Ny+1)
    vndmx   = vandermonde(Nx)
    vndmy   = vandermonde(Ny)
    patch   = Patch([1,1], vndmx'*coeffs*vndmy)
    ncoeffs = extractPatchCoeffs(patch)
    return L2norm(ncoeffs, coeffs, chebweights(Nx), chebweights(Ny))
end

function testinterpolatePatch(Nx::Int, Ny::Int)::Float64
    coeffs  = randn(Nx+1,Ny+1)
    vndmx   = vandermonde(Nx)
    vndmy   = vandermonde(Ny)
    patch   = Patch([1,1], vndmx'*coeffs*vndmy)
    npatch  = interpolatePatch(patch, Nx, Ny)
    return L2norm(npatch.value, patch.value, chebweights(Nx), chebweights(Ny))
end

function testinterpolatePatch(Nx::Int, Ny::Int, fn::Function)::Float64
    fpatch = Float64[fn(x,y) for x in chebgrid(Nx), y in chebgrid(Ny)]
    patch  = Patch([1,1], fpatch)
    npatch = interpolatePatch(patch, 2*Nx, 2*Ny)
    fpatch2N = Float64[fn(x,y) for x in chebgrid(2*Nx), y in chebgrid(2*Ny)]
    return L2norm(npatch.value, fpatch2N, chebweights(2*Nx), chebweights(2*Ny))
end

function testcalcPatch(bnd1fn::Function, bnd2fn::Function, rhsfn::Function, Nx::Int, Ny::Int)::Float64
    bndx = getPatchIC(bnd1fn, 0, Ny, 1, 1) 
    bndy = getPatchIC(bnd2fn, 1, Ny, 1, 1)
    rhs  = RHS(rhsfn, Nx, Ny)
    dop  = derivOP(Nx,Ny)
    bop  = boundaryOP(Nx,Ny) 
    npatch = calcPatch(bndx, bndy, rhs, dop, bop, [1,1])
    fpatch = Float64[bnd1fn(x) + bnd1fn(y) for x in chebgrid(Nx), y in chebgrid(Ny)]
    return L2norm(npatch.value, fpatch, chebgrid(Nx), chebgrid(Ny))
end

fpatch = Float64[(x-1)^2 + (y-1)^2 for x in 1:10, y in 1:15]

@test getPatchIC(x->x, 0, 4, 2, 1).value ≈ (chebgrid(4) + 1)/2
@test getPatchIC(x->x, 1, 6, 2, 1).value ≈ (chebgrid(6) + 1)/2
@test getPatchIC(x->x, 1, 9, 2, 2).value ≈ (chebgrid(9) - 1)/2
@test getPatchBnd(Patch([1,1], fpatch), 0).value == Float64[(x-1)^2 + 14^2 for x in 1:10]
@test getPatchBnd(Patch([1,1], fpatch), 1).value == Float64[(y-1)^2 + 9^2 for y in 1:15]
@test testextractPatchCoeffs(7,3) < 1e-14
@test testinterpolatePatch(12,5) < 1e-14
@test testinterpolatePatch(12, 13, (x,y)->x^2 + y^3 + y^2*x^3) < 1e-14
@test testinterpolatePatch(29, 31, (x,y)-> sin(pi*x) + exp(-y^2)) < 1e-14
@test testcalcPatch(x->(x^9-1), y->(y^9-1), (x,y)->0, 9, 9) < 1e-13 
