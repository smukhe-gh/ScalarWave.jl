#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Define operations for 1D spaces
#--------------------------------------------------------------------

import Base: maximum, minimum, range, length

# FIXME: maximum and minimum functions are swapped. That's the only thing that seems to work
maximum(S::Type{T}) where {T<:Space{Tag, N, max, min}} where {Tag, N, max, min}  = min
minimum(S::Type{T}) where {T<:Space{Tag, N, max, min}} where {Tag, N, max, min}  = max
order(S::Type{T})   where {T<:Space{Tag, N, max, min}} where {Tag, N, max, min}  = N 
length(S::Type{T})  where {T<:Space{Tag, N, max, min}} where {Tag, N, max, min}  = N + 1

function derivative(S::Type{T})::Operator{S} where {T<:GaussLobatto{Tag,N}} where {Tag, N}
    return Operator(S, [derivative(S, i, j) for i in 1:N, j in 1:N])
end

function integral(S::Type{T})::Operator{S} where {T<:GaussLobatto{Tag,N}} where {Tag,N}
    return Operator(S, Diagonal([integral(S, i) for i in 1:N+1]))
end

function identity(S::Type{T})::Operator{S} where {T<:GaussLobatto{Tag}} where {Tag}
    return Operator(S, UniformScaling(1))
end

function incomingboundary(S::Type{T})::Operator{S} where {T<:Cardinal{Tag, N}} where {Tag, N}
    return Operator(S, Diagonal([1, zeros(N)...])) 
end

function outgoingboundary(S::Type{T})::Operator{S} where {T<:Cardinal{Tag, N}} where {Tag, N}
    return Operator(S, Diagonal([zeros(N)..., 1])) 
end
