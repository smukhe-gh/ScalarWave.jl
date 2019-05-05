#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Define operations for 1D spaces
#--------------------------------------------------------------------

function Base. reshape(u::Field{S}) where {S}
    return reshape(u.value, (prod(size(S))))
end

function Base. reshape(A::Operator{S}) where {S}
    return reshape(A.value, (prod(size(S)), prod(size(S))))
end

function Base. reshape(::Type{S}, u::Array{T,1})::Field{S} where {T,S}
    return Field(S, reshape(u, size(S)))
end

function Base. reshape(::Type{S}, A::Array{T,2})::Operator{S} where {T,S}
    return Operator(S, reshape(A, (size(S)..., size(S)...)))
end

function Base. *(A::Operator{S}, u::Field{S})::Field{S} where {S}
    return reshape(S, reshape(A)*reshape(u))
end

function Base. *(A::Operator{S}, B::Operator{S})::Field{S} where {S}
    return reshape(S, reshape(A)*reshape(B))
end

function Base. *(u::Field{S}, B::Operator{S})::Field{S} where {S}
    return reshape(S, Diagonal(reshape(u))*reshape(B))
end

function Base. \(A::Operator{S}, u::Field{S})::Field{S} where {S}
    return reshape(S, reshape(A)\reshape(u))
end
