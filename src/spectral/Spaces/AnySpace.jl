#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Define operations for ND spaces
#--------------------------------------------------------------------

using LinearAlgebra

function Base. reshape(u::Field{S}) where {S}
    return reshape(u.value, (prod(size(u.space))))
end

function Base. reshape(A::Operator{S}) where {S}
    return reshape(A.value, (prod(size(A.space)), prod(size(A.space))))
end

function Base. reshape(space::S, u::Array{T,1})::Field{S} where {T,S}
    return Field(space, reshape(u, size(space)))
end

function Base. reshape(space::S, A::Array{T,2})::Operator{S} where {T,S}
    return Operator(space, reshape(A, (size(space)..., size(space)...)))
end

function Base. *(A::Operator{S}, u::Field{S})::Field{S} where {S}
    @assert range(A.space) == range(u.space)
    return reshape(A.space, reshape(A)*reshape(u))
end

function Base. *(A::Operator{S}, B::Operator{S})::Field{S} where {S}
    @assert range(A.space) == range(B.space)
    return reshape(A.space, reshape(A)*reshape(B))
end

function Base. *(u::Field{S}, B::Operator{S})::Field{S} where {S}
    @assert range(u.space) == range(B.space)
    return reshape(u.space, Diagonal(reshape(u))*reshape(B))
end

function Base. \(A::Operator{S}, u::Field{S})::Field{S} where {S}
    @assert range(A.space) == range(u.space)
    return reshape(A.space, reshape(A)\reshape(u))
end

function Base. isapprox(u::Field{S}, v::Field{S})::Bool where {S}
    @assert range(u.space) == range(v.space)
    return u.value ≈ v.value
end

