#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Define operations for 1D spaces
#--------------------------------------------------------------------

import Base: +, -, *, /, ==, ≈, typeof, range, length

# dimensions and shape
order(S::Type{T}) where {T<:Cardinal{Tag, N}} where {Tag, N}  = N
dim(S::Type{T}) where {T<:Cardinal{Tag, N}} where {Tag, N}    = 1
range(S::Type{T}) where {T<:Cardinal{Tag, N}} where {Tag, N}  = 1:N+1
length(S::Type{T}) where {T<:Cardinal{Tag, N}} where {Tag, N}    = N+1

# type of space
typeof(S::Type{T}) where {T<:GaussLobatto{Tag, N}} where {Tag, N} = Float64
typeof(S::Type{T}) where {T<:GaussLobatto{Tag, N}} where {Tag, N} = Rational{BigInt}

# field / field
+(A::Field{S}, B::Field{S}) where {S} = Field(S, A.value + B.value)
-(A::Field{S}, B::Field{S}) where {S} = Field(S, A.value - B.value)
*(A::Field{S}, B::Field{S}) where {S<:Cardinal{Tag,N}} where {Tag, N} = Field(S, A.value .* B.value)
/(A::Field{S}, B::Field{S}) where {S<:Cardinal{Tag,N}} where {Tag, N} = Field(S, A.value ./ B.value)
==(A::Field{S}, B::Field{S}) where {S} = (A.value == B.value)
≈(A::Field{S}, B::Field{S}) where {S} = (A.value ≈ B.value)

# operator / operator
+(A::Operator{S}, B::Operator{S}) where {S} = Operator(S, A.value + B.value)
-(A::Operator{S}, B::Operator{S}) where {S} = Operator(S, A.value - B.value)
==(A::Operator{S}, B::Operator{S}) where {S} = (A.value == B.value)
≈(A::Operator{S}, B::Operator{S}) where {S} = (A.value ≈ B.value)

function *(A::Operator{S}, B::Operator{S})::Operator{S} where {S}
    C = similar(A.value)
    for index in CartesianRange(size(C))
        C[index] = sum(A.value[index.I[1],k]*B.value[k,index.I[2]] for k in range(S))
    end
    return Operator(S, C)
end

# operator / field
function *(A::Operator{S}, u::Field{S})::Field{S} where {S}
    v = similar(u.value)
    for index in range(S)
        v[index] = sum(A.value[index,k]*u.value[k] for k in range(S))
    end
    return Field(S, v)
end

function *(u::Field{S}, A::Operator{S})::Operator{S} where {S}
    B = similar(A.value)
    for index in CartesianRange(size(B))
        B[index] = u.value[index.I[1]]*A.value[index.I[1], index.I[2]]
    end
    return Operator(S, B)
end

function solve(A::Operator{S}, u::Field{S})::Field{S} where {S}
    return Field(S, A.value \ u.value)
end

# evaluate a field at grid points
function Field(S::Type{Cardinal{Tag}}, umap::Function)::Field{S} where {Tag}
    value = zeros(typeof(S), length(S))
    for index in range(S)
        value[index] = umap(collocation(typeof(S), index, order(S)))
    end
    return Field(S, value)
end

# compute derivative
function derivative(S::Type{Cardinal{Tag}})::Operator{S} where {Tag}
    DS = zeros(typeof(S), length(S), length(S))
    for index in CartesianRange(size(DS))
        DS[index] = derivative(typeof(S), index.I[1], index.I[2], order(S))
    end
    return Operator(S, DS)
end

# compute identity matrix
function eye(S::Type{Cardinal{Tag}})::Operator{S} where {Tag}
    return Operator(S, eye(typeof(S), length(S)))
end

function boundary(S::Type{Cardinal{Tag}})::Operator{S} where {Tag}
    B = zeros(typeof(S), length(S))
    B[1] = B[end] = 1
    return Operator(S, diagm(vec(B)))
end

# map boundaries
function Boundary(S::Type{Cardinal{Tag}}, f::Function...)::Boundary{S} where {Tag},
    b      = zeros(typeof(S), length(S))
    b[1]   = f[1](collocation(typeof(S), 1, order(S)))
    b[end] = f[2](collocation(typeof(S), length(S), order(S)))
    typeof(bnd1val) <: typeof(S) ? b[1]   = bnd1val : error("Mapping doesn't preserve eltype. Aborting.")
    typeof(bnd2val) <: typeof(S) ? b[end] = bnd2val : error("Mapping doesn't preserve eltype. Aborting.")
    return Boundary(S, b)
end

# field / boundary
function +(u::Field{S}, b::Boundary{S})::Field{S} where {S}
    return Field(S, u.value + b.value)
end
