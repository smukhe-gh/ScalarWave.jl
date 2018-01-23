#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

mutable struct Patch
    loc::Array{Int,1}
    value::Array{Float64,2}
end

struct Boundary
    row::Array{Float64,1}
    col::Array{Float64,1}
end

function getPB(patch::Patch, s::Symbol)::Array{Float64,1}
    if s==:R
        boundary = patch.value[end, :]
    else
        boundary = patch.value[:, end]
    end
    return boundary
end

function setPB(patch::Patch, boundary::Boundary)::Patch
    patch.value[1, :] = boundary.row
    patch.value[:, 1] = boundary.col
    return patch
end

function calcPatch(loc::Array{Int,1}, BC::Boundary, operator::Array{Float64, 4})::Patch
    N = size(operator)[1] - 1
    patch = Patch(loc, zeros(N+1, N+1))
    setPB(patch, BC) 
    patch.value = shapeB(reshapA(operator)\reshapeB(patch.value)) 
    return patch
end
