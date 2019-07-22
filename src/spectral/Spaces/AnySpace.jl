#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Define operations for ND spaces
#--------------------------------------------------------------------

using LinearAlgebra
export ⊕, solve, cond, L1, L2

function Base. reshape(u::Field{S}) where {S}
    return reshape(u.value, (prod(size(u.space))))
end

function Base. reshape(A::Operator{S}) where {S}
    return reshape(A.value, (prod(size(A.space)), prod(size(A.space))))
end

function Base. reshape(space::S, u::AbstractArray{T,1})::Field{S} where {T,S}
    return Field(space, reshape(u, size(space)))
end

function Base. reshape(space::S, A::AbstractArray{T,2})::Operator{S} where {T,S}
    return Operator(space, reshape(A, (size(space)..., size(space)...)))
end

function Base. *(A::Operator{S}, u::Field{S})::Field{S} where {S}
    @assert range(A.space) == range(u.space)
    return reshape(A.space, reshape(A)*reshape(u))
end

function Base. *(A::Operator{S}, B::Operator{S})::Operator{S} where {S}
    @assert range(A.space) == range(B.space)
    return reshape(A.space, reshape(A)*reshape(B))
end

function Base. +(A::Operator{S}, B::Operator{S})::Operator{S} where {S}
    @assert range(A.space) == range(B.space)
    return reshape(A.space, reshape(A) + reshape(B))
end

function Base. -(A::Operator{S}, B::Operator{S})::Operator{S} where {S}
    @assert range(A.space) == range(B.space)
    return reshape(A.space, reshape(A) - reshape(B))
end

function Base. *(u::Field{S}, B::Operator{S})::Operator{S} where {S}
    @assert range(u.space) == range(B.space)
    return reshape(u.space, Diagonal(reshape(u))*reshape(B))
end

function Base. *(a::Number, B::Operator{S})::Operator{S} where {S}
    return reshape(B.space, eltype(reshape(B))(a).*reshape(B))
end

function solve(A::Operator{S}, u::Field{S})::Field{S} where {S}
    @assert range(A.space) == range(u.space)
    return reshape(A.space, reshape(A)\reshape(u))
end

function Base. isapprox(u::Field{S}, v::Field{S})::Bool where {S}
    @assert range(u.space) == range(v.space)
    return u.value ≈ v.value
end

function Base. *(u::Field{S}, v::Field{S})::Field{S} where {S}
    @assert range(u.space) == range(v.space)
    return Field(u.space, u.value.*v.value)
end
    
function Base. /(u::Field{S}, v::Field{S})::Field{S} where {S}
    @assert range(u.space) == range(v.space)
    return Field(u.space, u.value./v.value)
end

function  ⊕(A::Operator{S}, B::Operator{S})::Operator{S} where {S}
    @assert range(A.space) == range(B.space)
    ID = identity(A.space)
    return (ID - B)*A + B
end

function  ⊕(A::Operator{S}, B::Operator{S}, C::Operator{S})::Operator{S} where {S}
    @assert range(A.space) == range(B.space)
    ID = identity(A.space)
    return (ID - B)*A + B*C
end

function LinearAlgebra. cond(A::Operator{S}) where {S}
    return cond(reshape(A))
end

function LinearAlgebra. eigvals(A::Operator{S}) where {S}
    return sort(abs.(eigvals(reshape(A))))
end

function Base. +(u::Field{S}, v::Field{S})::Field{S} where {S}
    @assert range(u.space) == range(v.space)
    return reshape(u.space, reshape(u) + reshape(v))
end

function Base. -(u::Field{S}, v::Field{S})::Field{S} where {S}
    @assert range(u.space) == range(v.space)
    return reshape(u.space, reshape(u) - reshape(v))
end

function Base. ^(u::Field{S}, a::Number)::Field{S} where {S}
    return reshape(u.space, reshape(u).^a)
end

function Base. abs(u::Field{S})::Field{S} where {S}
    return Field(u.space, abs.(u.value))
end

function Base. maximum(u::Field{S})::Number where {S}
    return maximum(u.value)
end
    
function Base. *(a::Number, u::Field{S})::Field{S} where {S}
    return reshape(u.space, a.*reshape(u))
end

function Base. /(a::Number, u::Field{S})::Field{S} where {S}
    return reshape(u.space, a./reshape(u))
end

function Base. display(A::Operator{S}) where {S}
    display(reshape(A))
    println()
end

function L1(u::Field{S})::Number where {S}
    return maximum(abs(u))
end

function Base. sum(u::Field{S})::Number  where {S}
    return sum(u.value)
end

function L2(u::Field{S})::Number where {S}
    W = integral(u.space)
    return sqrt(sum(W*(u^2)))
end


