#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function pconv(bnd1::Function, bnd2::Function, Nx::Int, Ny::Int)::Float64
    rPatch = distribute(bnd1, bnd2, (x,y)->0, Nx, Ny, 1)[[1,1]]
    fPatch = Float64[bnd1(x) + bnd2(y) for x in chebgrid(4Nx, M, 1), y in chebgrid(4Ny, M, 1)]
    return L2norm(sPatch, interpolatePatch(rPatch, 4Nx, 4Ny).value, chebweights(4Nx), chebweights(4Ny))
end

function hconv(bnd1::Function, bnd2::Function, Nx::Int, Ny::Int, M::Int)::Float64
    dbase = distribute(bnd1, bnd2, (x,y)->0, Nx, Ny, M)
    sumL2 = 0.0
    for m in 1:M, n in 1:M
        sPatch = interpolatePatch(dbase[[m,n]], 2Nx, 2Ny).value
        fPatch = Float64[bnd1(x) + bnd2(y) for x in chebgrid(2Nx, M, m), y in chebgrid(2Ny, M, n)]
        sumL2 += (L2norm(fPatch, sPatch, chebweights(2Nx)/M, chebweights(2Ny)/M))^2
    end
    return sqrt(sumL2)
end

function convergence(bnd1::Function, bnd2::Function, spanN::Array{Int,1}, spanM::Arrat{Int,1})::Float64
    if size(spanM)[1] ==  1
        println("==> Starting p convergence")
        for p in spanN
            L2errornorm = pconv(bnd1, bnd2, p, p)
            println("p = $p np = $(spanM[1])  L2 error = $L2errornorm")
        end
    elseif size(spanN)[1] == 1
        println("==> Starting h convergence")
        L20 = 1
        for m in spanM
            L2errornorm = hconv(bnd1, bnd2, spanN[1], spanN[1], 2^m)
            println("p = $(spanN[1]), np = $(2^m)  L2 error = $L2errornorm ratio = $(L20/L2n)")
            L20 = L2n
        end
    end
end
