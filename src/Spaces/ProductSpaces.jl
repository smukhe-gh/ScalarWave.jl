#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Define operations and data structures for 1D spaces
#--------------------------------------------------------------------

import Base: +, -, *

struct ProductSpace{S1<:Space, S2<:Space} end

struct ProductSpaceOperator{S, D, T}
    space::Type{S}
    value::Array{T, D}
end

order{S1, S2}(PS::Type{ProductSpace{S1, S2}}) = (order(S1), order(S2)) 
dim{S1, S2}(PS::Type{ProductSpace{S1, S2}})   = dim(S1) + dim(S2)  
range{S1, S2}(PS::Type{ProductSpace{S1, S2}}) = CartesianRange((order(S1)+1, order(S2)+1))
size{S1, S2}(PS::Type{ProductSpace{S1, S2}})  = range(PS).stop.I

+{S}(A::ProductSpaceOperator{S}, B::ProductSpaceOperator{S}) = ProductSpaceOperator(S, A.value + B.value)
-{S}(A::ProductSpaceOperator{S}, B::ProductSpaceOperator{S}) = ProductSpaceOperator(S, A.value - B.value)

function Field{S1, S2}(PS::Type{ProductSpace{S1, S2}}, umap::Function)::Field{PS}
    value = zeros(size(PS))
    for index in range(PS) 
        value[index] = umap(chebx(index.I[1], order(S1)), chebx(index.I[2], order(S2))) 
    end
    return Field(PS, value)
end

function ⦼(A::Operator, B::Operator)::ProductSpaceOperator
    AB = zeros(size(A.space)..., size(B.space)..., 
               size(A.space)..., size(B.space)...)
    for index in CartesianRange(size(AB))
        i ,j ,k ,l = index.I
        AB[index]  = A.value[i,k]*B.value[j,l] 
    end
    return ProductSpaceOperator(ProductSpace{A.space, B.space}, AB)
end

function derivative{S1, S2}(PS::Type{ProductSpace{S1, S2}})     
    return (derivative(S1)⦼identity(S2), identity(S1)⦼derivative(S2))
end

function *{S}(A::ProductSpaceOperator{S}, B::ProductSpaceOperator{S})::ProductSpaceOperator{S} 
    C = similar(A.value)
    orders = order(S)
    for index in CartesianRange(size(C))
        i, j, k, l = index.I
        C[index]   = sum(A.value[i,j,m,n]*B.value[m,n,k,l] for m in 1:order(S)[1]+1, n in 1:order(S)[2]+1)
    end
    return ProductSpaceOperator(S, C)
end

function *{S}(A::ProductSpaceOperator{S}, B::Field{S})::Field{S} 
    C = similar(B.value)
    orders = order(S)
    for index in CartesianRange(size(C))
        i, j = index.I
        C[index]   = sum(A.value[i,j,m,n]*B.value[m,n] for m in 1:order(S)[1]+1, n in 1:order(S)[2]+1)
    end
    return Field(S, C)
end

function *{S}(A::Field{S}, B::ProductSpaceOperator{S})::ProductSpaceOperator{S} 
    C = similar(B.value)
    orders = order(S)
    for index in CartesianRange(size(C))
        i, j ,k ,l = index.I
        C[index]   = sum(B[i,j]*A.value[m,n,k,l] for m in 1:order(S)[1]+1, n in 1:order(S)[2]+1)
    end
    return Field(S, C)
end

function solve{S}(A::ProductSpaceOperator{S}, u::Field{S})::Field{S}
    return Field(S, reshapeL2H(reshapeH2L(A.value) \ reshapeH2L(u.value), size(S)))
end

function boundary{S1, S2}(PS::Type{ProductSpace{S1, S2}})     
    B = zeros(size(S1)..., size(S2)..., 
              size(S1)..., size(S2)...)
    for index in CartesianRange(size(B))
        k, i, j, l = index.I
        if  i==1 || k==1 || i == size(S1) || k == size(S2)
            B[index] = delta(i,j)*delta(k,l)
        end
    end
    return ProductSpaceOperator(ProductSpace{S1, S2}, B)
end

function Boundary{S1, S2}(PS::Type{ProductSpace{S1, S2}}, umap::Function)::Boundary{PS}
    value = zeros(size(PS))
    for index in range(PS) 
        i, j = index.I
        if  i==1 || j==1 || i == size(S1) || j == size(S2)
            value[index] = umap(chebx(index.I[1], order(S1)), chebx(index.I[2], order(S2))) 
        end
    end
    return Field(PS, value)
end
