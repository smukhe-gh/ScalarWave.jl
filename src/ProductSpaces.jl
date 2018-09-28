#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Define operations and data structures for 2D spaces
#--------------------------------------------------------------------

struct ProductSpace{S1<:Space, S2<:Space} end

struct Null end
struct Spatial end 

struct ProductSpaceOperator{S, D, T}
    space::Type{S}
    value::Array{T, D}
end

import Base: +, -, *, /, size, range, vec, sqrt, abs, zero, similar

order(PS::Type{ProductSpace{S1, S2}}) where {S1, S2} = (order(S2), order(S1)) 
dim(PS::Type{ProductSpace{S1, S2}}) where {S1, S2} = dim(S1) + dim(S2)  
range(PS::Type{ProductSpace{S1, S2}}) where {S1, S2} = CartesianRange((len(S2), len(S1)))
size(PS::Type{ProductSpace{S1, S2}}) where {S1, S2} = (len(S2), len(S1)) 

# FIXME: This is type-unstable? 
zero(u::Type{T}) where {T<:Field} = 0.0 
zero(::Type{ProductSpace{S1, S2}}) where {S1, S2}   = Field(ProductSpace{S1, S2}, (u,v)->0)  

function zero(::Type{Spatial}, PS::Type{ProductSpace{S1, S2}})::ProductSpaceOperator{ProductSpace{S1, S2}} where {S1, S2 <: GaussLobatto{Tag, N}}  where {Tag, N}  
    B = zeros(Float64, len(S2), len(S1))
    return ProductSpaceOperator(ProductSpace{S1, S2}, reshape(diagm(vec(B)), (len(S2), len(S1), len(S2), len(S1))))
end

#-----------------------------------------------------------------------------------
-(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, -A.value)
+(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, +A.value)
sqrt(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, sqrt.(A.value))
abs(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, abs.(A.value))

+(A::Field{ProductSpace{S1, S2}}, 
  B::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, A.value .+ B.value)
-(A::Field{ProductSpace{S1, S2}}, 
  B::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, A.value .- B.value)

*(a::T, 
  B::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag}, T<:Real} where {Tag} = Field(ProductSpace{S1, S2}, a .* B.value)
+(a::T, 
  B::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag}, T<:Real} where {Tag} = Field(ProductSpace{S1, S2}, a .+ B.value)
-(a::T, 
  B::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag}, T<:Real} where {Tag} = Field(ProductSpace{S1, S2}, a .- B.value)
/(a::T, 
  B::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag}, T<:Real} where {Tag} = Field(ProductSpace{S1, S2}, a ./ B.value)
/(B::Field{ProductSpace{S1, S2}}, a::T) where {S1, S2 <: Cardinal{Tag}, T<:Real} where {Tag} = Field(ProductSpace{S1, S2}, B.value ./ a)

*(A::Field{ProductSpace{S1, S2}}, 
  B::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, A.value .* B.value)
/(A::Field{ProductSpace{S1, S2}}, 
  B::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, A.value ./ B.value)

*(a::T, A::ProductSpaceOperator{PS}) where {T<:Real, PS} = ProductSpaceOperator(PS, a.*A.value)
/(A::ProductSpaceOperator{PS}, b::T) where {T<:Real, PS} = ProductSpaceOperator(PS, A.value ./b)

#-----------------------------------------------------------------------------------

+(A::ProductSpaceOperator{PS}, B::ProductSpaceOperator{PS}) where {PS} = ProductSpaceOperator(PS, A.value + B.value)
-(A::ProductSpaceOperator{PS}, B::ProductSpaceOperator{PS}) where {PS} = ProductSpaceOperator(PS, A.value - B.value)


function Field(PS::Type{ProductSpace{S1, S2}}, umap::Function)::Field{PS} where {S1, S2 <: GaussLobatto{Tag,N}} where {Tag, N}
    value = zeros(Float64, size(PS))
    for index in range(PS)
        value[index] = umap(collocation(Float64, index.I[1], order(S2)), 
                            collocation(Float64, index.I[2], order(S1))) 
    end
    return Field(PS, value)
end

function Field(PS::Type{ProductSpace{S1, S2}}, umap::Function)::Field{PS} where {S1, S2 <: Taylor{Tag,N}} where {Tag, N}
    value = zeros(Rational{BigInt}, size(PS))
    for index in range(PS) 
        value[index] = umap(collocation(Rational, index.I[1], order(S2)), 
                            collocation(Rational, index.I[2], order(S1))) 
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
    AB = zeros(Rational{BigInt}, len(S2), len(S1), len(S2), len(S1)) 
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
    ni = len(S2)
    nj = len(S1)
    n = ni * nj
    C = reshape(reshape(A.value, n, n) * reshape(B.value, n, n), ni, nj, ni, nj)
    return ProductSpaceOperator(ProductSpace{S1, S2}, C)
end

function *(A::ProductSpaceOperator{ProductSpace{S1,S2}}, 
           u::Field{ProductSpace{S1, S2}})::Field{ProductSpace{S1, S2}} where {S1, S2} 
    v = similar(u.value)
    for index in CartesianRange(size(v))
        i, ii = index.I
        v[index] = sum(A.value[i, ii, k, kk]*u.value[k, kk] for k in range(S2), kk in range(S1))
    end
    return Field(ProductSpace{S1, S2}, v)
end

function *(u::Field{PS}, A::ProductSpaceOperator{PS})::ProductSpaceOperator{PS} where {PS} 
    B = similar(A.value)
    for index in CartesianRange(size(B))
        i, j, ii, jj = index.I
        B[index]     = u.value[i, j]*A.value[index]
    end
    return ProductSpaceOperator(PS, B)
end

#--------------------------------------------------------------
# Spatial
#--------------------------------------------------------------
function boundary(::Type{Spatial}, PS::Type{ProductSpace{S1, S2}})::ProductSpaceOperator{ProductSpace{S1, S2}} where {S1, S2 <: GaussLobatto{Tag, N}}  where {Tag, N}  
    B = zeros(Float64, len(S2), len(S1))
    B[1, :] = B[:, 1] = B[end, :] = B[:, end] = 1
    return ProductSpaceOperator(ProductSpace{S1, S2}, reshape(diagm(vec(B)), (len(S2), len(S1), len(S2), len(S1))))
end

function boundary(::Type{Spatial}, PS::Type{ProductSpace{S1, S2}})::ProductSpaceOperator{ProductSpace{S1, S2}} where {S1, S2 <: Taylor{Tag, N}}  where {Tag, N}  
    B = zeros(Rational{BigInt}, len(S2), len(S1))
    B[1, :] = B[:, 1] = B[end, :] = B[:, end] = 1//1
    return ProductSpaceOperator(ProductSpace{S1, S2}, reshape(diagm(vec(B)), (len(S2), len(S1), len(S2), len(S1))))
end

#--------------------------------------------------------------
#Null 
#--------------------------------------------------------------
function boundary(::Type{Null}, PS::Type{ProductSpace{S1, S2}})::ProductSpaceOperator{ProductSpace{S1, S2}} where {S1, S2 <: GaussLobatto{Tag, N}}  where {Tag, N}  
    B = zeros(Float64, len(S2), len(S1))
    B[1, :] = B[:, 1] = 1
    return ProductSpaceOperator(ProductSpace{S1, S2}, reshape(diagm(vec(B)), (len(S2), len(S1), len(S2), len(S1))))

end

function boundary(::Type{Null}, PS::Type{ProductSpace{S1, S2}})::ProductSpaceOperator{ProductSpace{S1, S2}} where {S1, S2 <: Taylor{Tag, N}}  where {Tag, N}  
    B = zeros(Rational{BigInt}, len(S2), len(S1))
    B[1, :] = B[:, 1] =  1//1
    return ProductSpaceOperator(ProductSpace{S1, S2}, reshape(diagm(vec(B)), (len(S2), len(S1), len(S2), len(S1))))
end

function Boundary(PS::Type{ProductSpace{S1, S2}}, bmap::Function...)::Boundary{PS} where {S1, S2 <: GaussLobatto{Tag, N}}  where {Tag, N}  
    B = zeros(Float64, len(S2), len(S1))
    @assert length(bmap) == 4
    bnd1 = map(bmap[1], [collocation(Float64, i,  order(S2)) for i  in range(S2)])
    bnd2 = map(bmap[2], [collocation(Float64, ii, order(S1)) for ii in range(S1)])
    bnd3 = map(bmap[3], [collocation(Float64, j,  order(S2)) for j  in range(S2)])
    bnd4 = map(bmap[4], [collocation(Float64, jj, order(S1)) for jj in range(S1)])
    B[:,1] = bnd1 
    B[1,:] = bnd2 
    B[:, end] = bnd3 
    B[end, :] = bnd4 
    return Boundary(ProductSpace{S1, S2}, B)
end

function Boundary(PS::Type{ProductSpace{S1, S2}}, bmap::Function...)::Boundary{PS} where {S1, S2 <: Taylor{Tag, N}}  where {Tag, N}  
    B = zeros(Rational{BigInt}, len(S2), len(S1))
    bnd1 = map(bmap[1], [collocation(Rational, i,  order(S2)) for i  in range(S2)])
    bnd2 = map(bmap[2], [collocation(Rational, ii, order(S1)) for ii in range(S1)])
    bnd3 = map(bmap[3], [collocation(Rational, j,  order(S2)) for j  in range(S2)])
    bnd4 = map(bmap[4], [collocation(Rational, jj, order(S1)) for jj in range(S1)])
    eltype(bnd1) <: Rational ? B[:,1] = bnd1 : error("Mapping doesn't preserve eltype. Aborting") 
    eltype(bnd2) <: Rational ? B[1,:] = bnd2 : error("Mapping doesn't preserve eltype. Aborting") 
    eltype(bnd3) <: Rational ? B[:, end] = bnd3 : error("Mapping doesn't preserve eltype. Aborting") 
    eltype(bnd4) <: Rational ? B[end, :] = bnd4 : error("Mapping doesn't preserve eltype. Aborting") 
    return Boundary(ProductSpace{S1, S2}, B)
end

function +(u::Field{ProductSpace{S1,S2}}, b::Boundary{ProductSpace{S1,S2}})::Field{ProductSpace{S1,S2}} where {S1,S2}
    return Field(ProductSpace{S1,S2}, u.value + b.value)
end

function +(u::Real, b::Boundary{ProductSpace{S1,S2}})::Field{ProductSpace{S1,S2}} where {S1,S2}
    return Field(ProductSpace{S1,S2}, u + b.value)
end

function vec(u::Field{ProductSpace{S1, S2}})::Array{eltype(u.value),1} where {S1, S2} 
    return vec(u.value)
end

import Base.maximum

function maximum(u::Field{ProductSpace{S1, S2}}) where {S1, S2} 
    return maximum(u.value)
end

function vec(A::ProductSpaceOperator{ProductSpace{S1, S2}})::Array{eltype(A.value),2} where {S1, S2} 
    return reshape(A.value, (prod(size(ProductSpace{S1, S2})), prod(size(ProductSpace{S1, S2}))))
end

function shape(PS::Type{ProductSpace{S1, S2}}, u::Array{T,1})::Array{eltype(u),2} where {S1, S2, T}
    return reshape(u, size(PS))
end

function solve(A::ProductSpaceOperator{ProductSpace{S1, S2}}, u::Field{ProductSpace{S1, S2}})::Field{ProductSpace{S1, S2}} where {S1, S2}
    X =  vec(A) \ vec(u)
    return Field(ProductSpace{S1, S2}, shape(ProductSpace{S1, S2}, X))
end

import Base: eigvals, cond

function eigvals(A::ProductSpaceOperator{ProductSpace{S1, S2}}) where {S1, S2}
    return eigvals(vec(A))
end

function cond(A::ProductSpaceOperator{ProductSpace{S1, S2}}) where {S1, S2}
    return cond(vec(A))
end

function similar(u::Field{S, D, T})::Field{S,D,T} where {S, D, T}
    return Field(u.space, Array{T,D}(size(u.space)))
end

