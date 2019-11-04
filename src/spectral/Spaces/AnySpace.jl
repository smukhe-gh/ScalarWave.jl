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
    return reshape(u.space, reshape(u).*reshape(v))
end
    
function Base. /(u::Field{S}, v::Field{S})::Field{S} where {S}
    @assert range(u.space) == range(v.space)
    return reshape(u.space, reshape(u)./reshape(v))
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
    return reshape(u.space, abs.(reshape(u)))
end

function Base. maximum(u::Field{S})::Number where {S}
    return maximum(reshape(u))
end
    
function Base. minimum(u::Field{S})::Number where {S}
    return minimum(reshape(u))
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

function Base. display(u::Field{S}) where {S}
    display(u.value)
    println()
end

function L1(u::Field{S})::Number where {S}
    return maximum(abs(u))
end

function Base. sum(u::Field{S})::Number  where {S}
    return sum(reshape(u))
end

function L2(u::Field{S})::Number where {S}
    W = integral(u.space)
    return sqrt(sum(W*(u^2)))
end

function  ⊕(u::Field{S}, v::Field{S})::Field{S} where {S}
    @assert range(u.space) == range(v.space)
    I = identity(u.space)
    B = incomingboundary(u.space)
    return (I - B)*v + B*u
end

function Base. -(u::Field{S})::Field{S} where {S}
    return reshape(u.space, -reshape(u))
end

function Base. +(a::Number, u::Field{S})::Field{S} where {S}
    return reshape(u.space, a .+ reshape(u))
end

function Base. /(u::Field{S}, a::Number)::Field{S} where {S}
    return reshape(u.space, reshape(u)./a)
end

function Base. log(u::Field{S})::Field{S} where {S}
    return reshape(u.space, log.(reshape(u)))
end

function Base. log10(u::Field{S})::Field{S} where {S}
    return reshape(u.space, log10.(reshape(u)))
end

function Base. exp(u::Field{S})::Field{S} where {S}
    return reshape(u.space, exp.(reshape(u)))
end

function Base. sin(u::Field{S})::Field{S} where {S}
    return reshape(u.space, sin.(reshape(u)))
end

function Base. cos(u::Field{S})::Field{S} where {S}
    return reshape(u.space, cos.(reshape(u)))
end

function Base. sqrt(u::Field{S})::Field{S} where {S}
    return reshape(u.space, sqrt.(reshape(u)))
end

function Base. transpose(u::Field{S})::Field{S} where {S}
    return Field(u.space, transpose(u.value))
end
