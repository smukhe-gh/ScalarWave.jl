#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function testrestrictmodes1D(fn::Function, Nx::Int, M::Int)::Float64
    fpatch  = Float64[fn(x) for x in chebgrid(Nx)]
    vndmx   = vandermonde(Nx)
    frpatch = vndmx*restrictmodes!(inv(vndmx)*fpatch, M)
    w       = chebweights(Nx)
    L2error = sqrt(w'*(fpatch - frpatch).^2)
    return L2error
end

function testrestrictmodes2D(fn::Function, Nx::Int, Ny::Int, M::Int, N::Int)::Float64
    fpatch  = Float64[fn(x,y) for x in chebgrid(Nx), y in chebgrid(Ny)]
    vndmx   = vandermonde(Nx)
    vndmy   = vandermonde(Ny)
    frpatch = vndmx*restrictmodes!(inv(vndmx)*fpatch*inv(vndmy)', M, N)*vndmy'
    return L2norm(fpatch, frpatch, chebweights(Nx), chebweights(Ny))
end

function testprolongatemodes1D(fn::Function, Nx::Int, M::Int)::Float64
    fpatch  = Float64[fn(x) for x in chebgrid(Nx)]
    vndmx   = vandermonde(Nx)
    vndmMx  = pseudovandermonde(M, chebgrid(Nx))
    fppatch = vndmMx*prolongatemodes(inv(vndmx)*fpatch, M)
    w       = chebweights(Nx)
    return  sqrt(w'*(fpatch - fppatch).^2)
end

function testprolongatemodes2D(fn::Function, Nx::Int, Ny::Int, M::Int, N::Int)::Float64
    fpatch  = Float64[fn(x,y) for x in chebgrid(Nx), y in chebgrid(Ny)]
    vndmx   = vandermonde(Nx)
    vndmy   = vandermonde(Ny)
    vndmMx  = pseudovandermonde(M, chebgrid(Nx))
    vndmNy  = pseudovandermonde(N, chebgrid(Ny))
    fppatch = vndmMx*prolongatemodes(inv(vndmx)*fpatch*inv(vndmy)', M, N)*vndmNy'
    return L2norm(fpatch, fppatch, chebweights(Nx), chebweights(Ny))
end

function testprolongate1Dx(fn::Function, Nx::Int, M::Int)::Float64
    patch1D = Float64[fn(x) for x in chebgrid(Nx)]
    patchMD = prolongateOP(Nx, M)*patch1D
    wx      = chebweights(Nx)
    L2error = sum(sqrt(wx'*(patchMD[1+(m-1)*(Nx+1):m*(Nx+1)] - Float64[fn(x) for x in chebgrid(Nx, M, m)]).^2) for m in 1:M)
    return L2error 
end
 
function testprolongate1Dy(fn::Function, Ny::Int, M::Int)::Float64
    patch1D = Float64[fn(y) for y in chebgrid(Ny)]
    patchMD = patch1D'*prolongateOP(Ny, M)'
    wy      = chebweights(Ny)
    L2error = sum(sqrt(wy'*(patchMD'[1+(m-1)*(Ny+1):m*(Ny+1)] - Float64[fn(y) for y in chebgrid(Ny, M, m)]).^2) for m in 1:M)
    return L2error 
end

function testrestrict1Dx(fn::Function, Nx::Int, M::Int)::Float64
    patch1D = Float64[fn(x) for x in chebgrid(Nx)]
    patchMD = prolongateOP(Nx, M)*patch1D
    patchSD = restrictOP(Nx, M)*patchMD
    wx      = chebweights(Nx)
    L2error = sqrt(wx'*(patch1D - patchSD).^2)  
    return L2error 
end

function testrestrict1Dy(fn::Function, Ny::Int, M::Int)::Float64
    patch1D = Float64[fn(y) for y in chebgrid(Ny)]
    patchMD = patch1D'*prolongateOP(Ny, M)'
    patchSD = patchMD*restrictOP(Ny,M)'
    wy      = chebweights(Ny)
    L2error = sqrt(wy'*(patch1D - patchSD').^2)  
    return L2error 
end

function testrProlongateThenRestrictPatch(fn::Function, Nx::Int, Ny::Int, M::Int)::Float64
    sPatch   = Patch([1,1], Float64[fn(x,y) for x in chebgrid(Nx), y in chebgrid(Ny)])
    mPatch   = prolongatePatch(sPatch, M)
    m2sPatch = restrictPatch(mPatch)
    return L2norm(sPatch.value, m2sPatch.value, chebweights(Nx), chebweights(Ny))
end

@test testrestrictmodes1D(x->x^9, 12, 4) > testrestrictmodes1D(x->x^9, 12, 5)
@test testrestrictmodes2D((x,y)->x^9+y^8, 12, 12, 4, 5) > testrestrictmodes2D((x,y)->x^9+y^8, 12, 12, 7, 8)
@test testprolongatemodes1D(x->x^9, 10, 12) == testprolongatemodes1D(x->x^9, 10, 13)
@test testprolongatemodes2D((x,y)->x^2+y^2, 9, 9, 10, 10) == testprolongatemodes2D((x,y)->x^2+y^2, 9, 9, 12, 12) 
@test prolongateOP(4,2)*restrictOP(4,2)*prolongateOP(4,2) ≈ prolongateOP(4,2) 
@test restrictOP(4,2)*prolongateOP(4,2)*restrictOP(4,2) ≈ restrictOP(4,2)
@test (prolongateOP(4,2)*restrictOP(4,2))' ≈ prolongateOP(4,2)*restrictOP(4,2) 
@test (restrictOP(4,2)*prolongateOP(4,2))' ≈ restrictOP(4,2)*prolongateOP(4,2)
@test testprolongate1Dx(x->x.^3, 10, 3) < 1e-14
@test testprolongate1Dy(x->x.^5, 12, 2) < 1e-14
@test testrestrict1Dx(x->x.^5, 10, 4) < 1e-14
@test testrestrict1Dy(x->x.^2, 4, 2) < 1e-14
@test testrProlongateThenRestrictPatch((x,y)->x^2 + y^3, 12, 8, 2) < 1e-14
