#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2019
# Define operations for 3D spaces
#--------------------------------------------------------------------

import Base: maximum, minimum, range, size

maximum(S::Type{T}) where {T<:ProductSpace{S1, S2, S3}} where {S1, S2, S3} = (maximum(S1), maximum(S2), maximum(S3))
minimum(S::Type{T}) where {T<:ProductSpace{S1, S2, S3}} where {S1, S2, S3} = (minimum(S1), minimum(S2), minimum(S3))
order(S::Type{T})   where {T<:ProductSpace{S1, S2, S3}} where {S1, S2, S3} = (  order(S1),   order(S2),   order(S3))
size(S::Type{T})    where {T<:ProductSpace{S1, S2, S3}} where {S1, S2, S3} = (   size(S1),    size(S2),    size(S3)) 

function Base. kron(A::Operator{S1}, B::Operator{S2}, C::Operator{S3})::Operator{ProductSpace{S1, S2, S3}} where {S1, S2, S3}
    return reshape(ProductSpace{S1, S2, S3}, kron(kron(reshape(A), reshape(B)), reshape(C)))
end

function Base. identity(S::Type{T})::Operator{ProductSpace{S1, S2}} where {T<:ProductSpace{S1, S2, S3}} where {S1, S2, S3}
    return kron(identity(S1), identity(S2), identity(S3)) 
end

function incomingboundary(S::Type{T})::Operator{ProductSpace{S1,  S2, S3}} where {T<:ProductSpace{S1, S2, S3}} where {S1, S2, S3}
    return (kron(incomingboundary(S1), identity(S2), identity(S3)) + 
            kron(identity(S1), incomingboundary(S2), identity(S3)) + 
            kron(identity(S1), identity(S2), incomingboundary(S3))) 
end

function outgoingboundary(S::Type{T})::Operator{ProductSpace{S1,  S2}} where {T<:ProductSpace{S1, S2}} where {S1, S2}
    return (kron(outgoingboundary(S1), identity(S2), identity(S3)) + 
            kron(identity(S1), outgoingboundary(S2), identity(S3)) + 
            kron(identity(S1), identity(S2), outgoingboundary(S3))) 
end

function derivative(S::Type{T})::NTuple{2, Operator{ProductSpace{S1,  S2}}} where {T<:ProductSpace{S1, S2}} where {S1, S2}
    return (kron(derivative(S1), identity(S2), identity(S3)), 
            kron(identity(S1), derivative(S2), identity(S3)),
            kron(identity(S1), identity(S2), derivative(S3)))
end

function integral(S::Type{T})::Operator{ProductSpace{S1,  S2}} where {T<:ProductSpace{S1, S2}} where {S1, S2}
    return kron(integral(S1), integral(S2), integral(S3)) 
end
