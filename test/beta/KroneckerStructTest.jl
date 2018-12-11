#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 10-2018
# Test the structs for Kronecker products and the 
# the associated operations
#--------------------------------------------------------------------

import Base.getindex, Base. *, Base. +

abstract type AbstractOperator{S} end

struct Operator{S, T} <: AbstractOperator{S}
    space::S
    array::Array{T}
end

struct IdentityOperator{S} <: AbstractOperator{S} end

struct KroneckerProduct{S1, S2, S} <: AbstractOperator{S}
    O1::AbstractOperator{S1}
    O2::AbstractOperator{S2}
    function KroneckerProduct{S1, S2}() where {S1, S2}
        S = ProductSpace{S1, S2}
        new{S1, S2, S}()
    end
end

struct KroneckerAdd{S}
    O1::KroneckerProduct{S1, S2}
    O2::KroneckerProduct{S1, S2}
end

Base. +(x::AbstractOperator, y::AbstractOperator) = KroneckerAdd(x, y)

struct KroneckerMultiply{S1, S2}
    O1::KroneckerProduct{S1, S2}
    O2::KroneckerProduct{S1, S2}
end

Base. *(x::AbstractOperator, y::AbstractOperator) = KroneckerMultiply(x, y)

# compute the order of the space
order(S::Type{T}) where {T<:Cardinal{Tag, N}} where {Tag, N}  = N
function order(S::KroneckerProduct{S1, S2}) where {S1, S2}
    return (order(S1), order(S2))
end

# functions that create the structs
function kronecker(A::Union{Type{IdentityOperator{S1}}, Operator{S1}}, 
                   B::Union{Type{IdentityOperator{S2}}, Operator{S2}})::KroneckerProduct{S1, S2} where {S1, S2 <: Space}
    return KroneckerProduct(A, B)
end

function +(A::KroneckerProduct, B::KroneckerProduct)::KroneckerSum
    return KroneckerSum(A, B)
end

function *(A::KroneckerProduct{S1, S2}, B::KroneckerProduct{S1, S2})::KroneckerMultiply{S1, S2} where {S1, S2 <: Space}
    return KroneckerMultiply(A, B)
end

# getindex operations on the structs (evaulations are restricted to getindex calls)
function getindex(A::Operator{S}, indices...) where {S}
    i, j = indices
    return A.value[i,j]
end

function getindex(A::Type{IdentityOperator{S}}, indices...) where {S}
    i, j = indices
    @assert i <= order(S) + 1
    @assert j <= order(S) + 1
    if i==j
        return 1
    else 
        return 0
    end
end

function getindex(A::KroneckerProduct, indices...)
    i, j, k, l = indices
    return A.O1[i,k]*A.O2[k,l]
end

function getindex(A::KroneckerSum, indices...)
    i, j, k, l = indices
    return A.O1[i,j,k,l] + A.O2[i,j,k,l]
end

function getindex(A::KroneckerMultiply, indices...)
    i, j, k, l = indices
    return sum(A.O1[i,j,m,n]*A.O2[m,n,k,l] for m in order(A.O1)[2], n in order(A.O1)[1])
end

# operations between two constructs [inelegant]
function +(A::Union{KroneckerMultiply, KroneckerSum}, 
           B::Union{KroneckerMultiply, KroneckerSum})
    return KroneckerSum(A, B)
end

function *(A::Union{KroneckerMultiply, KroneckerSum}, 
           B::Union{KroneckerMultiply, KroneckerSum})
    return KroneckerProduct(A, B)
end

#--------------------------------------------------------------------
# Test
#--------------------------------------------------------------------

DU = derivative(GaussLobatto(U, 4))
DV = derivative(GaussLobatto(V, 4))
IU = IdentityOperator{GaussLobatto(U,4)}
IV = IdentityOperator{GaussLobatto(V,4)}

L2 = kronecker(DU, IV) + kronecker(IU, DV)
L1 = kronecker(DU, IV) * kronecker(IU, DV)
L  = L1 + L2

@show typeof(L1)
@show typeof(L2)
@show typeof(L)

