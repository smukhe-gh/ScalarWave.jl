#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function testdistribute(bnd1::Function, bnd2::Function, Nx::Int, Ny::Int, M::Int)::Float64
    dbase = distribute(bnd1, bnd2, (x,y)->0, Nx, Ny, M)
    sumL2 = 0.0
    for m in 1:M, n in 1:M
        sPatch = interpolatePatch(dbase[[m,n]], 2Nx, 2Ny).value
        fPatch = Float64[bnd1(x) + bnd2(y) for x in chebgrid(2Nx, M, m), y in chebgrid(2Ny, M, n)]
        sumL2 += (L2norm(fPatch, sPatch, chebweights(2Nx)/M, chebweights(2Ny)/M))^2
    end
    return sqrt(sumL2)
end

@test testdistribute(x->sin(pi*x), y->sin(pi*y), 24, 24, 1) < 1e-14
