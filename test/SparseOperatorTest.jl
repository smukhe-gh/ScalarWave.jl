#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 06-2019
# Test an implementation of Kronecker Structs a.k.a. sparse operators
# The algebra of sparse operators is a ring +, *
#--------------------------------------------------------------------
#    DenseOperator     DenseOperator
#         |                 |
#          -----------------
#                 | ⦼
#          SparseOperator{⦼}    SparseOperator{⦼}
#                 |                    |
#                  --------------------
#                           | +, *
#          SparseOperator{+}, SparseOperator{⦼}
#--------------------------------------------------------------------

abstract type Operator end

struct M end
struct ⦼ end

struct DenseOperator{S, D, T} <: Operator
    space::S
    value::AbstractArray{T, D}
end

struct SparseScalingOperator{S, T} <: Operator
    space::S
    value::AbstractArray{T,D}
end

struct SparseOperator{S, T} <: Operator
    O1::Operator
    O2::Operator
end

function Base. *(A::DenseOperator{S}, B::DenseOperator{S})::DenseOperator{S} where {S}
    return DenseOperator(A.space, A.value*B.value)
end

function Base. *(A::SparseOperator{S, ⦼}, B::SparseOperator{S, ⦼})::SparseOperator{S, ⦼} where {S}
    return SparseOperator{S, ⦼}(A.O1*B.O1, A.O2*B.O2)
end

function Base. *(u::Field{S}, B::SparseOperator{S, ⦼})::SparseOperator{S, ⦼} where {S}
    return SparseOperator{S, *}(SparseScalingOperator(u.space, u.value), B)
end

function Base. collect(A::SparseOperator{S, ⦼})::AbstractArray  where {S}
    @warn "Colleting a Kronecker struct is an expensive operation. Avoid this"
    return kron(A.O1.value, A.O2.value)
end

function Base. +(A::SparseOperator{S, ⦼}, B::SparseOperator{S, ⦼})::SparseOperator{S, +} where {S}
    return SparseOperator{S, +}(A, B)
end

function Base. +(A::SparseOperator{S, +}, B::SparseOperator{S, ⦼})::SparseOperator{S, +} where {S}
    return SparseOperator{S, +}(A, B)
end

function Base. getindex(A::DenseOperator{S}, i::T, j::T)::Number where {S} where {T<:Int}
    return A.value[i,j]
end

function Base. getindex(A::SparseScalingOperator{S}, i::T, j::T)::Number where {S} where {T<:Int}
    return A.value[i,k]
end

function Base. getindex(A::SparseOperator{S, ⦼}, i::T, j::T, k::T, l::T)::Number where {S} where {T<:Int}
    return A.O1[i, k]*A.O2[j, l]
end

function Base. getindex(A::SparseOperator{S, +}, i::T, j::T, k::T, l::T)::Number where {S} where {T<:Int}
    return A.O1[i,j,k,l] + A.O2[i,j,k,l]
end

function Base. getindex(A::SparseOperator{S, *}, i::T, j::T, k::T, l::T)::Number where {S} where {T<:Int}
    return A.O1[i,k] + A.O2[i,j,k,l]
end

function Base. getindex(A::SparseOperator{S, ⦼}, i::T, j::T)::Number where {S} where {T<:Int}
    (m, n) = size(A.O1.value)    
    (k, l) = size(A.O2.value)    
    return A.O1[cld(i, k), cld(j, l)] * A.O2[(i - 1) % k + 1, (j - 1) % l + 1]
end

function Base. getindex(A::SparseOperator{S, +}, i::T, j::T)::Number where {S} where {T<:Int}
    return A.O1[i,j] + A.O2[i,j]
end

function Base. getindex(A::SparseOperator{S, *}, i::T, j::T)::Number where {S} where {T<:Int}
    return A.O1[i,i]*A.O2[i,j]
end

#--------------------------------------------------------------------
# Now test this.
#--------------------------------------------------------------------

function test2index(A::SparseOperator{S, ⦼})::AbstractArray  where {S}
    @warn "Colleting a Kronecker struct is an expensive operation. Avoid this"
    (m, n) = size(A.O1.value)    
    (p, q) = size(A.O2.value)    
    value  = zeros(m*p, n*q)
    for index in CartesianIndices(value)
        value[index] = A[index.I...]
    end
    return reshape(value, (m*p, n*q))
end

function test4index(A::SparseOperator{S, ⦼})::AbstractArray  where {S}
    @warn "Colleting a Kronecker struct is an expensive operation. Avoid this"
    (m, n) = size(A.O1.value)    
    (p, q) = size(A.O2.value)    
    value  = zeros(m, p, n, q)
    for index in CartesianIndices(value)
        value[index] = A[index.I...]
    end
    return reshape(value, (m*p, n*q))
end

Sx = ChebyshevGL{M, 2, BigFloat}(-1, 1)
Sy = ChebyshevGL{M, 2, BigFloat}(-3, 5)
Dx = DenseOperator(Sx, derivative(Sx).value)
Dy = DenseOperator(Sy, derivative(Sy).value)

D2x = Dx*Dx 
D2y = Dy*Dy

D1 = SparseOperator{M, ⦼}(Dx, Dy) 
D2 = D1 + D1 + D1
D3 = D1 * D1

@show typeof(D1)
@show typeof(D2)
@show typeof(D3)

@test D3.O1.value ≈ D2x.value
@test D3.O2.value ≈ D2y.value
@test D2[1, 1, 2, 1] ≈ 3*D1[1, 1, 2, 1]
@test collect(D1) ≈ kron(Dx.value, Dy.value)
