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

function testprojectonPatchBnd(fn::Function, Nx::Int)::Float64
    pfbnd = projectonPatchBnd(fn, Nx) 
    efbnd = Float64[fn(x) for x in chebgrid(Nx)]
    w     = chebweights(Nx)  
    L2error = sqrt(w'*(pfbnd-efbnd).^2/w'*(efbnd).^2)
    return L2error 
end

function testprojectonPatch(fn::Function, Nx::Int, Ny::Int)::Float64
    pfpatch = projectonPatch(fn, Nx, Ny)
    efpatch = Float64[fn(x,y) for x in chebgrid(Nx), y in chebgrid(Ny)]
    return L2norm(pfpatch, efpatch, chebweights(Nx), chebweights(Ny))
end

fpatch = Float64[(x-1)^2 + (y-1)^2 for x in 1:10, y in 1:15]

@test getPatchIC(4, 3, 2, [1,1], x->x, 0).value == (chebgrid(4) + 1)/2
@test getPatchIC(5, 6, 2, [2,1], x->x, 1).value == (chebgrid(6) + 1)/2
@test getPatchIC(4, 9, 2, [1,2], x->x, 1).value == (chebgrid(9) - 1)/2
@test getPatchBnd(Patch([1,1], fpatch), 0).value == Float64[(x-1)^2 + 14^2 for x in 1:10]
@test getPatchBnd(Patch([1,1], fpatch), 1).value == Float64[(y-1)^2 + 9^2 for y in 1:15]
@test testextractPatchCoeffs(7,3) < 1e-14
@test testinterpolatePatch(12,5) < 1e-14
@test testinterpolatePatch(8, 9, (x,y)->x^2 + y^3 + y^2*x^3) < 1e-14
@test testinterpolatePatch(19, 21, (x,y)-> sin(pi*x) + exp(-y^2)) < 1e-14
@test_broken testprojectonPatchBnd(x->x^2+1, 4) < 1e-14
@test_broken testprojectonPatch((x,y)->x^2+y^3+x^3*y^2, 8, 9) < 1e-14



