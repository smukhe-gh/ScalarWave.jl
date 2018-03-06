#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function testdistribute(Nx::Int, Ny::Int, M::Int)::Float64
    dbase = distribute(x->sin(pi*x), y->sin(pi*y), (x,y)->0, Nx, Ny, M)
    sumL2 = 0.0
    for m in 1:M, n in 1:M
        sPatch = interpolatePatch(dbase[[m,n]], 2Nx, 2Ny).value
        fPatch = Float64[sin(pi*x) + sin(pi*y) for x in chebgrid(2Nx, M, m), y in chebgrid(2Ny, M, n)]
        sumL2 += L2norm(fPatch, sPatch, chebweights(2Nx)/M, chebweights(2Ny)/M)
    end
    return sumL2
end

@test testdistribute(24, 24, 1) < 1e-14
@test_broken testdistribute(22, 22, 2) < 1e-14
