#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 09-2019
# Utilities for interfacing with NLsolve solvers
#--------------------------------------------------------------------

export reshapeFromTuple, reshapeToTuple
export reshapeFromTuple2E, reshapeToTuple2E
export symmetrize

function reshapeFromTuple(U::NTuple{3, Field})
    return vcat(reshape(U[1]), reshape(U[2]), reshape(U[3]))
end

function reshapeToTuple(space::S, x::Array{T,1})::NTuple{3, Field}  where {S, T}
    U = reshape(x, :, 3)
    return (reshape(space, U[:, 1]), reshape(space, U[:, 2]), reshape(space, U[:, 3]))
end

function reshapeFromTuple2E(U::NTuple{2, Field})
    return vcat(reshape(U[1]), reshape(U[2]))
end

function reshapeToTuple2E(space::S, x::Array{T,1})::NTuple{2, Field}  where {S, T}
    U = reshape(x, :, 2)
    return (reshape(space, U[:, 1]), reshape(space, U[:, 2]))
end

function symmetrize(u::Field{S})::Field{S} where {S}
    return (u + transpose(u))/2 
end


