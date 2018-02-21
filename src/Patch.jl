#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function getPatchBnd(patch::Patch, s::Symbol)::Boundary
    if s==:R
        boundary = Boundary(:R, patch.value[end, :])
    elseif s==:C
        boundary = Boundary(:C, patch.value[:, end])
    else 
        error("Unknown symbol passed")
    end
    return boundary
end

function calcPatch(loc::Array{Int,1}, bndx::Boundary, bndy::Boundary, operator::Array{Float64, 4})::Patch
    Nx = size(operator)[1] - 1
    Ny = size(operator)[2] - 1
    B  = zeros(Nx+1, Ny+1)

    if bnd0.value[1] != bnd1.value[1]
        error("Inconsistent boundary conditions.")
    else
        B[1, :] = bndx.value
        B[:, 1] = bndy.value
    end
    
    return Patch(loc, shapeL2H(shapeH2L(operator) \ shapeH2L(B))) 
end

function extractPatchCoeffs(patch::Patch)::Array{Float64,2}
    fnodal = patch.value
    N      = size(fnodal)[1] - 1
    invndm = inv(vandermonde(N,chebgrid(N)))
    fmodal = zeros(N+1, N+1)
    for i in 1:N+1, j in 1:N+1
        elem = 0.0
        for m in 1:N+1, n in 1:N+1
            elem = elem + invndm[i,m]*fnodal[m,n]*invndm[j,n]
        end  
        fmodal[i,j] = elem
    end
    return fmodal
end
