#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Define operations for 1D spaces
#--------------------------------------------------------------------

import Base: +, -, *, /, ==, ≈, range, length, minimum, maximum, real, imag

# field / field
+(A::Field{S}, B::Field{S}) where {S} = Field(S, A.value + B.value)
-(A::Field{S}, B::Field{S}) where {S} = Field(S, A.value - B.value)
*(A::Field{S}, B::Field{S}) where {S<:Cardinal{Tag,N}} where {Tag, N} = Field(S, A.value .* B.value)
/(A::Field{S}, B::Field{S}) where {S<:Cardinal{Tag,N}} where {Tag, N} = Field(S, A.value ./ B.value)
==(A::Field{S}, B::Field{S}) where {S} = (A.value == B.value)
≈(A::Field{S}, B::Field{S}) where {S} = (A.value ≈ B.value)

# scalar / field
+(a::T, B::Field{S}) where {T<:Real, S} = Field(S, a .+ B.value)
-(a::T, B::Field{S}) where {T<:Real, S} = Field(S, a .- B.value)
*(a::T, B::Field{S}) where {T<:Real, S<:Cardinal{Tag,N}} where {Tag, N} = Field(S, a .* B.value)
/(a::T, B::Field{S}) where {T<:Real, S<:Cardinal{Tag,N}} where {Tag, N} = Field(S, a ./ B.value)

# operator / operator
+(A::Operator{S}, B::Operator{S}) where {S} = Operator(S, A.value + B.value)
-(A::Operator{S}, B::Operator{S}) where {S} = Operator(S, A.value - B.value)
==(A::Operator{S}, B::Operator{S}) where {S} = (A.value == B.value)
≈(A::Operator{S}, B::Operator{S}) where {S} = (A.value ≈ B.value)

function *(A::Operator{S}, B::Operator{S})::Operator{S} where {S}
    C = similar(A.value)
    for index in CartesianIndices(size(C))
        C[index] = sum(A.value[index.I[1],k]*B.value[k,index.I[2]] for k in range(S))
    end
    return Operator(S, C)
end

# operator / field
function *(u::Field{S}, A::Operator{S})::Operator{S} where {S}
    B = similar(A.value)
    for index in CartesianIndices(size(B))
        B[index] = u.value[index.I[1]]*A.value[index.I[1], index.I[2]]
    end
    return Operator(S, B)
end

function solve(A::Operator{S}, u::Field{S})::Field{S} where {S}
    return Field(S, A.value \ u.value)
end

# evaluate a field at grid points
function Field(S::Type{T}, umap::Function)::Field{S} where {T<:Cardinal{Tag, N}} where {Tag, N}
    value = zeros(size(S))
    for index in CartesianIndices(value) 
        value[index] = umap(collocation(S, index.I[1]))
    end
    return Field(S, value)
end

# create a zero coefficent field
function Field(S::Type{T})::Field{S} where {T<:Space} where {Tag, N}
    value = zeros(length(S))
    return Field(S, value)
end

# field / boundary
function +(u::Field{S}, b::Boundary{S})::Field{S} where {S}
    return Field(S, u.value + b.value)
end

# Define multiplication rules to extract operators from J
function Base. *(A::Operator{S}, u::Symbol)::Operator{S} where {S}
    return A
end

function Base. *(v::Field{S}, u::Symbol)::Operator{S} where {S}
    return v*eye(S)
end

function Base. *(A::Operator{S}, u::Int)::Operator{S} where {S}
    @assert u == 0 "Right mutiplication with an Int is only allowed to represent \n action on a zero vector field"
    return (A-A)                 
end

function Base. *(v::Field{S}, u::Int)::Operator{S} where {S}
    @assert u == 0 "Right mutiplication with an Int is only allowed to represent \n action on a zero vector field"
    return (eye(S) - eye(S))                
end

function ⊙(A::Operator{S}, B::Operator{S})::Operator{S} where {S}
    I = eye(S)
    return (I-B)*A + B*B
end

function ⊙(u::Field{S}, b::Boundary{S})::Field{S} where {S}
    I = eye(S)
    return (I-B)*u + B*b
end
