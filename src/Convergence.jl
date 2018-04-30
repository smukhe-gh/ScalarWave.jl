#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function sconv(bnd1::Function, bnd2::Function, Nx::Int, Ny::Int, M::Int)::Float64
    dbase = fdistribute(bnd1, bnd2, (x,y)->0, Nx, Ny, M)
    sumL2 = 0.0
    for m in 1:M, n in 1:M
        sPatch = interpolatePatch(dbase[[m,n]], 4Nx, 4Ny).value
        fPatch = Float64[bnd1(x) + bnd2(y) for x in chebgrid(4Nx, M, m), y in chebgrid(4Ny, M, n)]
        sumL2 += (L2norm(fPatch, sPatch, chebweights(4Nx)/M, chebweights(4Ny)/M))^2
    end
    return sqrt(sumL2)
end

function showconv(bnd1::Function, bnd2::Function, maxN::Int, maxM::Int, h::Int)
    if maxM > 1
        L2 = 1
        println("---------------------------------------------------------------------------------")
        println("==> h-convergence")
        println("---------------------------------------------------------------------------------")
        for m in 0:maxM
            L2error = sconv(bnd1, bnd2, maxN, maxN, h^m)
            m == 0 ? R = 0 : R =  L2/L2error
            @printf("p = %e | np = %e | R = %e | L2 error = %e \n", maxN, h^m, R, L2error)
            L2 = L2error
        end
    else
        println("---------------------------------------------------------------------------------")
        println("==> p-convergence")
        println("---------------------------------------------------------------------------------")
        L2 = 1
        for n in 1:maxN
            L2error = sconv(bnd1, bnd2, n, n, maxM)
            n == 1 ? R = 0 : R =  L2/L2error
            @printf("p = %e | np = %e | R = %e | L2 error = %e \n", n, maxM, R, L2error)
            L2 = L2error
        end
    end
end
