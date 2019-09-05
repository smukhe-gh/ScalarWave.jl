#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2019
# Define operations for 1D spaces
#--------------------------------------------------------------------

using LinearAlgebra
import Base: maximum, minimum, range, size, ndims
export order, identity, incomingboundary, outgoingboundary

ndims(space::S)   where {S<:Cardinal{Tag, N, T}} where {Tag, N, T} = 1 
minimum(space::S) where {S<:Cardinal{Tag, N, T}} where {Tag, N, T} = space.min
maximum(space::S) where {S<:Cardinal{Tag, N, T}} where {Tag, N, T} = space.max
order(space::S)   where {S<:Cardinal{Tag, N, T}} where {Tag, N, T} = N-1 
size(space::S)    where {S<:Cardinal{Tag, N, T}} where {Tag, N, T} = N 
range(space::S)   where {S<:Cardinal{Tag, N, T}} where {Tag, N, T} = (minimum(space), maximum(space))

function Base. identity(space::S)::Operator{S} where {S<:Cardinal{Tag, N, T}} where {Tag, N, T}
    return Operator(space, Diagonal(ones(T, size(space))))
end

function outgoingboundary(space::S)::Operator{S} where {S<:Cardinal{Tag, N, T}} where {Tag, N, T}
    return Operator(space, Diagonal([1, zeros(T, order(space))...])) 
end

function incomingboundary(space::S)::Operator{S} where {S<:Cardinal{Tag, N, T}} where {Tag, N, T}
    return Operator(space, Diagonal([zeros(T, order(space))..., 1])) 
end

function derivative(space::S)::Operator{S} where {S<:Cardinal{Tag, N, T}} where {Tag, N, T}
    return Operator(space, [derivative(space, i, j) for i in 1:size(space), j in 1:size(space)])
end

function integral(space::S)::Operator{S} where {S<:Cardinal{Tag, N, T}} where {Tag, N, T}
    return Operator(space, Diagonal([integral(space, i) for i in 1:size(space)]))
end

function Field(space::S, umap::Function)::Field{S} where {S<:Cardinal{Tag, N, T}} where {Tag, N, T}
    value = zeros(T, size(space))
    for index in CartesianIndices(value)
        value[index] = umap(collocation(space, index.I[1]))
    end
    return Field(space, value)
end
