#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 09-2019
# Utilities for nonlinear solver
#--------------------------------------------------------------------

export reshapeFromTuple, reshapeToTuple 

function reshapeFromTuple(U::NTuple{3, Field})
    return vcat(reshape(U[1]), reshape(U[2]), reshape(U[3]))
end

function reshapeToTuple(space::S, x::Array{T,1})::NTuple{3, Field}  where {S, T}
    U = reshape(x, :, 3)
    return (reshape(space, U[:, 1]), reshape(space, U[:, 2]), reshape(space, U[:, 3]))
end

