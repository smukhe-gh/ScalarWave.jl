#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2019
# Define operations for 2D spaces
#--------------------------------------------------------------------

import Base: maximum, minimum, range, size, ndims

ndims(PS::T)   where {T<:ProductSpace{S1, S2}} where {S1, S2} = 2 
minimum(PS::T) where {T<:ProductSpace{S1, S2}} where {S1, S2} = (minimum(PS.S1), minimum(PS.S2))
maximum(PS::T) where {T<:ProductSpace{S1, S2}} where {S1, S2} = (maximum(PS.S1), maximum(PS.S2))
order(PS::T)   where {T<:ProductSpace{S1, S2}} where {S1, S2} = (order(PS.S1), order(PS.S2)) 
size(PS::T)    where {T<:ProductSpace{S1, S2}} where {S1, S2} = (size(PS.S1),   size(PS.S2)) 
range(PS::T)   where {T<:ProductSpace{S1, S2}} where {S1, S2} = (minimum(PS),   maximum(PS)) 

function Base. kron(A::Operator{S1}, B::Operator{S2})::Operator{ProductSpace{S2, S1}} where {S2, S1}
    return reshape(ProductSpace(B.space, A.space), kron(reshape(A), reshape(B)))
end

function Base. identity(PS::ProductSpace{S1, S2})::Operator{ProductSpace{S1, S2}} where {S1, S2}
    return kron(identity(PS.S2), identity(PS.S1)) 
end

function incomingboundary(PS::ProductSpace{S1, S2})::Operator{ProductSpace{S1,  S2}} where {S1, S2}
    return kron(identity(PS.S2), incomingboundary(PS.S1)) ⊕ kron(incomingboundary(PS.S2), identity(PS.S1))  
end

function outgoingboundary(PS::ProductSpace{S1, S2})::Operator{ProductSpace{S1,  S2}} where {S1, S2}
    return kron(identity(PS.S2), outgoingboundary(PS.S1)) ⊕ kron(outgoingboundary(PS.S2), identity(PS.S1))  
end

function derivative(PS::ProductSpace{S1, S2})::NTuple{2, Operator{ProductSpace{S1,  S2}}} where {S1, S2}
    return (kron(identity(PS.S2), derivative(PS.S1)), kron(derivative(PS.S2), identity(PS.S1)))
end

function integral(PS::ProductSpace{S1, S2})::Operator{ProductSpace{S1,  S2}}  where {S1, S2}
    return kron(integral(PS.S2), integral(PS.S1)) 
end

function Field(PS::ProductSpace{S1, S2}, umap::Function)::Field{ProductSpace{S1, S2}} where {S1, S2 <: Cardinal{Tag, N, T}} where {Tag, N, T}
    value = zeros(T, size(PS))
    for index in CartesianIndices(value)
        value[index] = umap(collocation(PS.S1, index.I[1]),
                            collocation(PS.S2, index.I[2]))
    end
    return Field(PS, value)
end
