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

# XXX: Convert these to element-wise expressions?
function extractPatchCoeffs{T<:Array{Float64,1}}(patch::Patch, xvec::T, yvec::T)::Array{Float64,2}
    N      = size(patch.value)[1] - 1
    fnodal = patch.value
    fmodal = inv(vandermonde(N,xvec))*fnodal*inv(vandermonde(N,xvec)')
    return fmodal
end

# XXX: Convert these to element-wise expressions?
function interpolatePatch{T<:Array{Float64,1}}(patch::Patch, xvec::T, yvec::T,  xinterp::T, yinterp::T)::Patch
    N      = size(patch.value)[1] - 1
    fmodal = extractPatchCoeffs(patch, xvec, yvec)
    fnodal = vandermonde(N,xinterp)*fmodal*vandermonde(N,yinterp)' 
    return Patch(patch.loc, fnodal)
end
