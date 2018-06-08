#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function sconv(bnd1::Function, bnd2::Function, frhs::Function, asol::Function, 
                                                                   Nx::Int, Ny::Int, M::Int)::Float64
    dbase = fdistribute(bnd1, bnd2, frhs, Nx, Ny, M)
    sumL2 = 0.0
    for m in 1:M, n in 1:M
        sPatch = interpolatePatch(dbase[[m,n]], 4Nx, 4Ny).value
        fPatch = Float64[asol(x,y) for x in chebgrid(4Nx, M, m), y in chebgrid(4Ny, M, n)]
        sumL2 += (L2norm(fPatch, sPatch, chebweights(4Nx)/M, chebweights(4Ny)/M))^2
    end
    return sqrt(sumL2)
end

function showconv(bnd1::Function, bnd2::Function, frhs::Function, analyticsol::Function, 
                                                               maxmodes::Int, maxlevels::Int, h::Int)
    if maxlevels > 1
        L2 = 1
        p  = maxmodes
        println("----------------------------------------------------------------------------------")
        println("==> h-convergence")
        println("----------------------------------------------------------------------------------")
        for level in 0:maxlevels
            L2error = sconv(bnd1, bnd2, frhs, analyticsol, p, p, h^level)
            level == 0 ? R = 0 : R =  L2/L2error
            @printf("p = %s | np = %s | h^(p+1) = %s | R = %e | L2 error = %e \n", 
                    lpad(p,4," "), lpad(h^level,4," "),  lpad(h^(p+1),4," "), R, L2error)
            L2 = L2error
        end
    else
        println("----------------------------------------------------------------------------------")
        println("==> p-convergence")
        println("----------------------------------------------------------------------------------")
        L2 = 1
        h  = 1
        level = 0
        for p in 1:maxmodes
            L2error = sconv(bnd1, bnd2, frhs, analyticsol, p, p, h^level)
            p == 1 ? R = 0 : R =  L2/L2error
            @printf("p = %s | np = %s | h^(p+1) = %s | R = %e | L2 error = %e \n", 
                    lpad(p,4," "), lpad(h^level,4," "),  lpad(h^(p+1),4," "), R, L2error)
            L2 = L2error
        end
    end
end
