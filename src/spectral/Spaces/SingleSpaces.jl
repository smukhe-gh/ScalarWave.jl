#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Define operations for 1D spaces
#--------------------------------------------------------------------

import Base: +, -, *, /, ==, ≈, range, length, minimum, maximum, real, imag

# Helper functions for creating spaces
GaussLobatto(Tag, N, max, min) = GaussLobatto{Tag ,N, max, min} 
Chebyshev(Tag, N, max, min)    = Chebyshev{Tag ,N, max, min}
Taylor(Tag, N, max, min)       = Taylor{Tag ,N, max, min}

# Default space is 1 to -1
GaussLobatto(Tag, N)           = GaussLobatto{Tag ,N, 1, -1} 
Chebyshev(Tag, N)              = Chebyshev{Tag ,N, 1, -1}
Taylor(Tag, N)                 = Taylor{Tag ,N, 1, -1}

# Real and Imaginary part of a field
real(A::Field{S}) where {S} = Field(S, real(A.value))
imag(A::Field{S}) where {S} = Field(S, imag(A.value))

# Pull out the Tag
Tag(S::Type{T}) where {T<:GaussLobatto{Tag, N, max, min}} where {Tag, N, max, min}  = Tag

# Pull out the minimum and the maximum coordinate bounds
minimum(S::Type{T}) where {T<:GaussLobatto{Tag, N, max, min}} where {Tag, N, max, min}  = min
maximum(S::Type{T}) where {T<:GaussLobatto{Tag, N, max, min}} where {Tag, N, max, min}  = max

minimum(S::Type{T}) where {T<:Chebyshev{Tag, N, max, min}} where {Tag, N, max, min}  = min
maximum(S::Type{T}) where {T<:Chebyshev{Tag, N, max, min}} where {Tag, N, max, min}  = max

minimum(S::Type{T}) where {T<:Taylor{Tag, N, max, min}} where {Tag, N, max, min}  = min
maximum(S::Type{T}) where {T<:Taylor{Tag, N, max, min}} where {Tag, N, max, min}  = max

# dimensions and shape
order(S::Type{T}) where {T<:Cardinal{Tag, N}} where {Tag, N}  = N
dim(S::Type{T}) where {T<:Cardinal{Tag, N}} where {Tag, N}    = 1
range(S::Type{T}) where {T<:Cardinal{Tag, N}} where {Tag, N}  = 1:N+1
length(S::Type{T}) where {T<:Cardinal{Tag, N}} where {Tag, N} = N+1
length(S::Type{T}) where {T<:Galerkin{Tag, N}} where {Tag, N} = N+1

# for Galerkin
order(S::Type{T}) where {T<:Galerkin{Tag, N}} where {Tag, N}  = N

# type of space
spacetype(S::Type{T}) where {T<:GaussLobatto{Tag, N}} where {Tag, N} = Float64
spacetype(S::Type{T}) where {T<:Taylor{Tag, N}} where {Tag, N} = Rational{BigInt}

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
function *(A::Operator{S}, u::Field{S})::Field{S} where {S}
    v = similar(u.value)
    for index in range(S)
        v[index] = sum(A.value[index,k]*u.value[k] for k in range(S))
    end
    return Field(S, v)
end

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
    value = zeros(spacetype(S), length(S))
    for index in range(S)
        value[index] = umap(collocation(S, index))
    end
    return Field(S, value)
end

# create a zero coefficent field
function Field(S::Type{T})::Field{S} where {T<:Space} where {Tag, N}
    value = zeros(length(S))
    return Field(S, value)
end

# compute derivative
function derivative(S::Type{T})::Operator{S} where {T<:Cardinal{Tag, N}} where {Tag, N}
    DS = zeros(spacetype(S), length(S), length(S))
    for index in CartesianIndices(size(DS))
        DS[index] = derivative(S, index.I[1], index.I[2])
    end
    return Operator(S, DS)
end

# compute integral and the only operation defined on it. 
function integral(S::Type{T})::IntegrationOperator{S} where {T<:Cardinal{Tag, N}} where {Tag, N}
    W = diagm(0=>[integral(S, i) for i in range(S)])
    return IntegrationOperator(S, W)
end

function *(W::IntegrationOperator{S}, u::Field{S})::Real where {S}
    return sum(W.value*u.value) 
end

# compute identity matrix
function eye(S::Type{T})::Operator{S} where {T<:Cardinal{Tag, N}} where {Tag, N}
    return Operator(S, Matrix{spacetype(S)}(I, length(S), length(S)))
end

function boundary(S::Type{T})::Operator{S} where {T<:Cardinal{Tag, N}} where {Tag, N}
    B = zeros(spacetype(S), length(S))
    B[1] = B[end] = 1
    return Operator(S, diagm(0 => vec(B)))
end

# map boundaries
function Boundary(S::Type{T}, f::Function...)::Boundary{S} where {T<:Cardinal{Tag, N}} where {Tag, N}
    b      = zeros(spacetype(S), length(S))
    bnd1val = f[1](collocation(S, 1))
    bnd2val = f[2](collocation(S, length(S)))
    typeof(bnd1val) <: spacetype(S) ? b[1]   = bnd1val : error("Mapping doesn't preserve eltype. Aborting.")
    typeof(bnd2val) <: spacetype(S) ? b[end] = bnd2val : error("Mapping doesn't preserve eltype. Aborting.")
    return Boundary(S, b)
end

# field / boundary
function +(u::Field{S}, b::Boundary{S})::Field{S} where {S}
    return Field(S, u.value + b.value)
end

# shape and reshape operators
function Base. vec(u::Field{S})::Array{eltype(u.value),1} where {S<:GaussLobatto{Tag,N}} where {Tag, N}
    return u.value
end

function Base. vec(A::Operator{S})::Array{eltype(A.value),2} where {S<:GaussLobatto{Tag,N}} where {Tag, N}
    return A.value
end
  
function shape(S::Type{T}, u::Array{Float64,1})::Field{S} where {T<:GaussLobatto{Tag,N}} where {Tag,N}
    return Field(S, u)
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
