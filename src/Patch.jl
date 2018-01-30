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
    fmodal = inv(vandermonde(N,x))*fnodal*inv(vandermonde(N,x)')
    return fmodal
end

function interpolatePatch(patch::Patch, x::Array{Float64,1}, y::Array{Float64,1})::Patch
    N      = size(patch.value)[1] - 1
    fmodal = extractPatchCoeffs(patch)
    fnodal = vandermonde(N,x)'*fmodal*vandermonde(N,y) 
    return Patch(patch.loc, fnodal)
end
