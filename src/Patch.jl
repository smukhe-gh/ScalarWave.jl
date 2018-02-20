#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

struct Patch
    loc::Array{Int,1}
    value::Array{Float64,2}
end

struct Boundary
    kind::Symbol
    value::Array{Float64,1}
end

function getPB(patch::Patch, s::Symbol)::Boundary
    if s==:R
        boundary = Boundary(:R, patch.value[end, :])
    elseif s==:C
        boundary = Boundary(:C, patch.value[:, end])
    else 
        error("Unknown symbol passed")
    end
    return boundary
end

function calcPatch(loc::Array{Int,1}, bnd0::Boundary, bnd1::Boundary, operator::Array{Float64, 4})::Patch
    N = size(operator)[1] - 1
    B = zeros(N+1, N+1)
    if bnd0.value[1] != bnd1.value[1]
        error("Inconsistent boundary conditions.")
    else
        B[1, :] = bnd0.value
        B[:, 1] = bnd1.value
    end
    return Patch(loc, shapeB(reshapeA(operator) \ reshapeB(B))) 
end

function extractPatchCoeffs(patch::Patch)::Array{Float64,2}
    fnodal = patch.value
    N      = size(fnodal)[1] - 1
    x      = Float64[chebx(i, N) for i in 1:N+1]
    invndm = inv(vandermonde(N,x))    
    fmodal = zeros(N+1, N+1)
    # TODO: Test loop
    for m in 1:N+1, j in 1:N+1
        elem = 0.0
        for k in 1:N+1, j in 1:N+1
            elem = elem + invndm[m,j]*fmodal[j,k]*invndm[n,k]
        end  
        fmodal[m,n] = elem
    end
    #fmodal = inv(vandermonde(N,x))*fnodal*inv(vandermonde(N,x)')
    return fmodal
end
