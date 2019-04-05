#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Define operations for 1D spaces
#--------------------------------------------------------------------

import Base: maximum, minimum, range, size

# FIXME: maximum and minimum functions are swapped. 
maximum(S::Type{T}) where {T<:Space{Tag, N, max, min}} where {Tag, N, max, min} = min
minimum(S::Type{T}) where {T<:Space{Tag, N, max, min}} where {Tag, N, max, min} = max
order(S::Type{T})   where {T<:Space{Tag, N, max, min}} where {Tag, N, max, min} = N 
size(S::Type{T})    where {T<:Space{Tag, N, max, min}} where {Tag, N, max, min} = N + 1 

function Base. identity(S::Type{T})::Operator{S} where {T<:Space{Tag}} where {Tag}
    return Operator(S, Diagonal(ones(size(S))))
end

function incomingboundary(S::Type{T})::Operator{S} where {T<:Space{Tag}} where {Tag}
    return Operator(S, Diagonal([1, zeros(order(S))...])) 
end

function outgoingboundary(S::Type{T})::Operator{S} where {T<:Space{Tag}} where {Tag}
    return Operator(S, Diagonal([zeros(order(S))..., 1])) 
end

function derivative(S::Type{T})::Operator{S} where {T<:GaussLobatto{Tag,N}} where {Tag, N}
    return Operator(S, [derivative(S, i, j) for i in 1:size(S), j in 1:size(S)])
end

function integral(S::Type{T})::Operator{S} where {T<:GaussLobatto{Tag,N}} where {Tag,N}
    return Operator(S, Diagonal([integral(S, i) for i in 1:size(S)]))
end

