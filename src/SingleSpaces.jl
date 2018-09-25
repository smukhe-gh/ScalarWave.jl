#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Define operations and data structures for 1D spaces
#--------------------------------------------------------------------

struct GaussLobatto{Tag ,N} <: Cardinal{Tag, N} end
struct Taylor{Tag ,N} <: Cardinal{Tag, N} end
struct Chebyshev{Tag ,N} <: Galerkin{Tag, N} end

struct Field{S, D, T}
    space::Type{S}
    value::Array{T, D}
end

struct Boundary{S, D, T}
    space::Type{S}
    value::Array{T, D}
end

struct Operator{S, D, T}
    space::Type{S}
    value::Array{T, D}
end

import Base: range, identity 
import Base: +, -, *, /, ==, ≈

order(S::Type{T}) where {T<:Cardinal{Tag, N}} where {Tag, N}  = N 
dim(S::Type{T}) where {T<:Cardinal{Tag, N}} where {Tag, N}    = 1 
range(S::Type{T}) where {T<:Cardinal{Tag, N}} where {Tag, N}  = 1:N+1 
len(S::Type{T}) where {T<:Cardinal{Tag, N}} where {Tag, N}    = N+1 

+(A::Field{S}, B::Field{S}) where {S} = Field(S, A.value + B.value)
-(A::Field{S}, B::Field{S}) where {S} = Field(S, A.value - B.value)
*(A::Field{S}, B::Field{S}) where {S<:Cardinal{Tag,N}} where {Tag, N} = Field(S, A.value .* B.value)
/(A::Field{S}, B::Field{S}) where {S<:Cardinal{Tag,N}} where {Tag, N} = Field(S, A.value ./ B.value)

+(A::Operator{S}, B::Operator{S}) where {S} = Operator(S, A.value + B.value)
-(A::Operator{S}, B::Operator{S}) where {S} = Operator(S, A.value - B.value)

==(A::Field{S}, B::Field{S}) where {S} = (A.value == B.value)
==(A::Operator{S}, B::Operator{S}) where {S} = (A.value == B.value)
≈(A::Field{S}, B::Field{S}) where {S} = (A.value ≈ B.value)
≈(A::Operator{S}, B::Operator{S}) where {S} = (A.value ≈ B.value)

function Field(S::Type{GaussLobatto{Tag,N}}, umap::Function)::Field{S} where {Tag, N}
    value = zeros(Float64, len(S))
    for index in range(S) 
        value[index] = umap(collocation(Float64, index, N)) 
    end
    return Field(S, value)
end

function Field(S::Type{Taylor{Tag,N}}, umap::Function)::Field{S} where {Tag, N}
    value = zeros(Rational{BigInt}, len(S))
    for index in range(S) 
        value[index] = umap(collocation(Rational, index, N)) 
        typeof(value[index]) <: Rational ? 1 : error("The mapping doens't preserve eltype. Aborting.")
    end
    return Field(S, value)
end

function derivative(S::Type{GaussLobatto{Tag,N}})::Operator{S} where {Tag, N}
    DS = zeros(len(S), len(S)) 
    for index in CartesianRange(size(DS))
        DS[index] = derivative(Float64, index.I[1], index.I[2], N) 
    end
    return Operator(S, DS)
end

function derivative(S::Type{Taylor{Tag,N}})::Operator{S} where {Tag, N}
    DS = zeros(Rational{BigInt}, len(S), len(S)) 
    for index in CartesianRange(size(DS))
        DS[index] = derivative(Rational, index.I[1], index.I[2], N) 
    end
    return Operator(S, DS)
end

function identity(S::Type{GaussLobatto{Tag,N}})::Operator{S} where {Tag, N}
    return Operator(S, eye(Float64, len(S)))
end

function identity(S::Type{Taylor{Tag,N}})::Operator{S} where {Tag, N}
    return Operator(S, eye(Rational{BigInt}, len(S)))
end

function boundary(S::Type{GaussLobatto{Tag,N}})::Operator{S} where {Tag, N}
    B = zeros(Float64, len(S)) 
    B[1] = B[end] = 1
    return Operator(S, diagm(vec(B)))
end

function boundary(S::Type{Taylor{Tag,N}})::Operator{S} where {Tag, N}
    B = zeros(Rational{BigInt}, len(S)) 
    B[1] = B[end] = 1//1
    return Operator(S, diagm(vec(B)))
end

function *(A::Operator{S}, B::Operator{S})::Operator{S} where {S} 
    C = similar(A.value)
    for index in CartesianRange(size(C))
        C[index] = sum(A.value[index.I[1],k]*B.value[k,index.I[2]] for k in range(S))
    end
    return Operator(S, C)
end

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

function Boundary(S::Type{GaussLobatto{Tag,N}}, f::Function...)::Boundary{S} where {Tag, N}
    b      = zeros(Float64, len(S))
    b[1]   = f[1](collocation(Rational, 1, order(S))) 
    b[end] = f[2](collocation(Rational, len(S), order(S))) 
    return Boundary(S, b) 
end

function Boundary(S::Type{Taylor{Tag,N}}, f::Function...)::Boundary{S} where {Tag, N}
    b      = zeros(Rational{BigInt}, len(S))
    bnd1val = f[1](collocation(Rational, 1, order(S))) 
    bnd2val = f[2](collocation(Rational, len(S), order(S))) 
    typeof(bnd1val) <: Rational ? b[1]   = bnd1val : error("Mapping doesn't preserve eltype. Aborting.")
    typeof(bnd2val) <: Rational ? b[end] = bnd2val : error("Mapping doesn't preserve eltype. Aborting.")
    return Boundary(S, b) 
end

function +(u::Field{S}, b::Boundary{S})::Field{S} where {S}
    return Field(S, u.value + b.value)
end

function solve(A::Operator{S}, u::Field{S})::Field{S} where {S}
    return Field(S, A.value \ u.value)
end
