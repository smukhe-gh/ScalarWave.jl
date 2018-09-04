#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Define operations and data structures for 2D spaces
#--------------------------------------------------------------------

struct ProductSpace{S1<:Space, S2<:Space} end

struct ProductSpaceOperator{S, D, T}
    space::Type{S}
    value::Array{T, D}
end

import Base: +, -, *, size, range

order(PS::Type{ProductSpace{S1, S2}}) where {S1, S2} = (order(S1), order(S2)) 
dim(PS::Type{ProductSpace{S1, S2}}) where {S1, S2} = dim(S1) + dim(S2)  
range(PS::Type{ProductSpace{S1, S2}}) where {S1, S2} = CartesianRange((len(S1), len(S2)))
size(PS::Type{ProductSpace{S1, S2}}) where {S1, S2} = (len(S1), len(S2)) 

+(A::ProductSpaceOperator{S}, B::ProductSpaceOperator{S}) where {S} = ProductSpaceOperator(S, A.value + B.value)
-(A::ProductSpaceOperator{S}, B::ProductSpaceOperator{S}) where {S} = ProductSpaceOperator(S, A.value - B.value)

function Field(PS::Type{ProductSpace{S1, S2}}, umap::Function)::Field{PS} where {S1, S2 <: GaussLobatto{Tag,N}} where {Tag, N}
    value = zeros(Rational, size(PS))
    for index in range(PS) 
        value[index] = umap(collocation(Float64, index.I[1], order(S1)), 
                            collocation(Float64, index.I[2], order(S2))) 
    end
    return Field(PS, value)
end

function Field(PS::Type{ProductSpace{S1, S2}}, umap::Function)::Field{PS} where {S1, S2 <: Taylor{Tag,N}} where {Tag, N}
    value = zeros(Rational, size(PS))
    for index in range(PS) 
        value[index] = umap(collocation(Rational, index.I[1], order(S1)), 
                            collocation(Rational, index.I[2], order(S2))) 
    end
    return Field(PS, value)
end

function ⦼(A::Operator{S1}, B::Operator{S2})::ProductSpaceOperator{ProductSpace{S1, S2}} where {S1, S2 <: GaussLobatto{Tag, N}} where {Tag, N}
    AB = zeros(Float64, len(S2), len(S1), len(S2), len(S1)) 
    for index in CartesianRange(size(AB))
        i ,j ,ii ,jj = index.I
        AB[index]    = A.value[j,jj]*B.value[i,ii] 
    end
    return ProductSpaceOperator(ProductSpace{S1, S2}, AB)
end

function ⦼(A::Operator{S1}, B::Operator{S2})::ProductSpaceOperator{ProductSpace{S1, S2}} where {S1, S2 <: Taylor{Tag, N}} where {Tag, N}
    AB = zeros(Rational, len(S2), len(S1), len(S2), len(S1)) 
    for index in CartesianRange(size(AB))
        i ,j ,ii ,jj = index.I
        AB[index]    = A.value[j,jj]*B.value[i,ii] 
    end
    return ProductSpaceOperator(ProductSpace{S1, S2}, AB)
end

function derivative(PS::Type{ProductSpace{S1, S2}}) where {S1, S2}     
    return (derivative(S1) ⦼ identity(S2), identity(S1) ⦼ derivative(S2))
end 

function *(A::ProductSpaceOperator{ProductSpace{S1,S2}}, 
           B::ProductSpaceOperator{ProductSpace{S1,S2}})::ProductSpaceOperator{ProductSpace{S1, S2}} where {S1, S2}
    C = similar(A.value)
    for index in CartesianRange(size(C))
        i, j, ii, jj = index.I
        C[index]     = sum(A.value[i,j,k,kk]*B.value[k,kk,ii,jj] for k in range(S2), kk in range(S1))    
    end
    return ProductSpaceOperator(ProductSpace{S1, S2}, C)
end

function *(A::ProductSpaceOperator{ProductSpace{S1,S2}}, 
           u::Field{ProductSpace{S1, S2}})::Field{ProductSpace{S1, S2}} where {S1, S2} 
    v = similar(B.value)
    for index in CartesianRange(size(C))
        i, ii = index.I
        C[index] = sum(A.value[i, ii, k, kk]*u.value[k, kk] for k in range(S2), kk in range(S1))
    end
    return Field(ProductSpace{S1, S2}, v)
end

function *(u::Field{PS}, A::ProductSpaceOperator{PS})::ProductSpaceOperator{PS} where {PS} 
    B = similar(B.value)
    for index in CartesianRange(size(C))
        i, j, ii, jj = index.I
        C[index]     = B[i, ii]*A.value[i, ii, j, jj]
    end
    return Field(S, C)
end

function boundary(PS::Type{ProductSpace{S1, S2}}) where {S1, S2 <: GaussLobatto{Tag, N}}  where {Tag, N}  
    B = zeros(Float64, len(S1), len(S2))
    B[1, :] = B[:, 1] = B[end, :] = B[:, end] = 1
    return ProductSpaceOperator(ProductSpace{S1, S2}, diagm(vec(B)))
end

function boundary(PS::Type{ProductSpace{S1, S2}}) where {S1, S2 <: Taylor{Tag, N}}  where {Tag, N}  
    B = zeros(Rational, len(S1), len(S2))
    B[1, :] = B[:, 1] = B[end, :] = B[:, end] = 1//1
    return ProductSpaceOperator(ProductSpace{S1, S2}, diagm(vec(B)))
end

function Boundary(PS::Type{ProductSpace{S1, S2}}, bmap::Function...)::Boundary{PS} where {S1, S2 <: GaussLobatto{Tag, N}}  where {Tag, N}  
    B = zeros(Float64, len(S1), len(S2))
    B[1, :]   = [bmap[1](collocation(Rational, ii, order(S1)) for ii in order(S1))]
    B[:, 1]   = [bmap[2](collocation(Rational, i,  order(S2)) for i  in order(S2))]
    B[:, end] = [bmap[3](collocation(Rational, jj, order(S1)) for jj in order(S1))]
    B[end, :] = [bmap[4](collocation(Rational, j,  order(S2)) for j  in order(S2))]
    return ProductSpaceOperator(ProductSpace{S1, S2}, B)
end

function Boundary(PS::Type{ProductSpace{S1, S2}}, bmap::Function...)::Boundary{PS} where {S1, S2 <: Taylor{Tag, N}}  where {Tag, N}  
    B = zeros(Rational, len(S1), len(S2))
    bnd1 = [bmap[1](collocation(Rational, ii, order(S1)) for ii in order(S1))]
    bnd2 = [bmap[2](collocation(Rational, i,  order(S2)) for i  in order(S2))]
    bnd3 = [bmap[3](collocation(Rational, jj, order(S1)) for jj in order(S1))]
    bnd4 = [bmap[4](collocation(Rational, j,  order(S2)) for j  in order(S2))]
    eltype(bnd1) <: Rational ? B[1,:] = bnd1 : error("Mapping doesn't preserve eltype. Aborting") 
    eltype(bnd2) <: Rational ? B[:,1] = bnd2 : error("Mapping doesn't preserve eltype. Aborting") 
    eltype(bnd3) <: Rational ? B[end,:]  = bnd3 : error("Mapping doesn't preserve eltype. Aborting") 
    eltype(bnd4) <: Rational ? B[:, end] = bnd4 : error("Mapping doesn't preserve eltype. Aborting") 
    return ProductSpaceOperator(ProductSpace{S1, S2}, B)
end

function solve(A::ProductSpaceOperator{PS}, u::Field{PS})::Field{PS} where {PS}
    return Field(S, reshapeL2H(reshapeH2L(A.value) \ reshapeH2L(u.value), size(PS)))
end
