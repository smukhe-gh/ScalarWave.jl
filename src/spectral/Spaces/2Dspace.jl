#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2019
# Define operations for 2D spaces
#--------------------------------------------------------------------

import Base: maximum, minimum, range, size

maximum(S::Type{T}) where {T<:ProductSpace{S1, S2}} where {S1, S2} = (maximum(S1), maximum(S2))
minimum(S::Type{T}) where {T<:ProductSpace{S1, S2}} where {S1, S2} = (minimum(S1), minimum(S2))
order(S::Type{T})   where {T<:ProductSpace{S1, S2}} where {S1, S2} = (  order(S1),   order(S2))
size(S::Type{T})    where {T<:ProductSpace{S1, S2}} where {S1, S2} = (   size(S1),    size(S2)) 

function Base. kron(A::Operator{S1}, B::Operator{S2})::Operator{ProductSpace{S1, S2}} where {S1, S2}
    return reshape(ProductSpace{S1, S2}, kron(reshape(A), reshape(B)))
end

function Base. identity(S::Type{T})::Operator{ProductSpace{S1, S2}} where {T<:ProductSpace{S1, S2}} where {S1, S2}
    return kron(identity(S1), identity(S2)) 
end

function incomingboundary(S::Type{T})::Operator{ProductSpace{S1,  S2}} where {T<:ProductSpace{S1, S2}} where {S1, S2}
    return kron(identity(S1), incomingboundary(S2)) + kron(incomingboundary(S1), identity(S2))  
end

function outgoingboundary(S::Type{T})::Operator{ProductSpace{S1,  S2}} where {T<:ProductSpace{S1, S2}} where {S1, S2}
    return kron(identity(S1), outgoingboundary(S2)) + kron(outgoingboundary(S1), identity(S2))  
end

function derivative(S::Type{T})::NTuple{2, Operator{ProductSpace{S1,  S2}}} where {T<:ProductSpace{S1, S2}} where {S1, S2}
    return (kron(derivative(S1), identity(S2)), kron(identity(S1), derivative(S2)))
end

function integral(S::Type{T})::Operator{ProductSpace{S1,  S2}} where {T<:ProductSpace{S1, S2}} where {S1, S2}
    return kron(integral(S1), integral(S2)) 
end
