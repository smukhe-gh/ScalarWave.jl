#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Define basic math and array operations for 2D spaces
#--------------------------------------------------------------------

import Base: +, -, *, /, size, range, vec, zero, similar, maximum, minimum

# identity element
zero(u::Type{T}) where {T<:Field} = 0.0
zero(::Type{ProductSpace{S1, S2}}) where {S1, S2} = Field(ProductSpace{S1, S2}, (u,v)->0)
zero(::Type{Null}, ::Type{ProductSpace{S1, S2}}) where {S1, S2} = 0*boundary(Null, ProductSpace{S1, S2})

# allocate memory for a similar field
similar(u::Field{S, D, T}) where {S, D, T} = Field(u.space, Array{T,D}(undef, size(u.space)))

# unary operators
-(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag, N}} where {Tag,N} = Field(ProductSpace{S1, S2}, -A.value)
+(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag, N}} where {Tag,N} = Field(ProductSpace{S1, S2}, +A.value)

# field / field
+(A::Field{ProductSpace{S1, S2}},
  B::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag, N}} where {Tag,N} = Field(ProductSpace{S1, S2}, A.value .+ B.value)
-(A::Field{ProductSpace{S1, S2}},
  B::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag, N}} where {Tag,N} = Field(ProductSpace{S1, S2}, A.value .- B.value)
*(A::Field{ProductSpace{S1, S2}},
  B::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag, N}} where {Tag,N} = Field(ProductSpace{S1, S2}, A.value .* B.value)
/(A::Field{ProductSpace{S1, S2}},
  B::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag, N}} where {Tag,N} = Field(ProductSpace{S1, S2}, A.value ./ B.value)

# scalar / field
+(a::T,
  B::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag, N}, T<:Real} where {Tag,N} = Field(ProductSpace{S1, S2}, a .+ B.value)
-(a::T,
  B::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag, N}, T<:Real} where {Tag,N} = Field(ProductSpace{S1, S2}, a .- B.value)
*(a::T,
  B::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag, N}, T<:Real} where {Tag,N} = Field(ProductSpace{S1, S2}, a .* B.value)
/(a::T,
  B::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag, N}, T<:Real} where {Tag,N} = Field(ProductSpace{S1, S2}, a ./ B.value)
/(B::Field{ProductSpace{S1, S2}}, a::T) where {S1, S2 <: Cardinal{Tag, N}, T<:Real} where {Tag,N} = Field(ProductSpace{S1, S2}, B.value ./ a)

+(a::T,
  B::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Space{Tag}, T<:Complex} where {Tag} = Field(ProductSpace{S1, S2}, a .+ B.value)
-(a::T,
  B::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Space{Tag}, T<:Complex} where {Tag} = Field(ProductSpace{S1, S2}, a .- B.value)
*(a::T,
  B::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Space{Tag}, T<:Complex} where {Tag} = Field(ProductSpace{S1, S2}, a .* B.value)

# scalar / operator
*(a::T, A::ProductSpaceOperator{PS}) where {T<:Real, PS} = ProductSpaceOperator(PS, a.*A.value)
/(A::ProductSpaceOperator{PS}, b::T) where {T<:Real, PS} = ProductSpaceOperator(PS, A.value ./b)

# operator / operator
+(A::ProductSpaceOperator{PS}, B::ProductSpaceOperator{PS}) where {PS} = ProductSpaceOperator(PS, A.value + B.value)
-(A::ProductSpaceOperator{PS}, B::ProductSpaceOperator{PS}) where {PS} = ProductSpaceOperator(PS, A.value - B.value)

function *(A::ProductSpaceOperator{ProductSpace{S1,S2}},
           B::ProductSpaceOperator{ProductSpace{S1,S2}})::ProductSpaceOperator{ProductSpace{S1, S2}} where {S1, S2}
    ni = length(S2)
    nj = length(S1)
    n = ni * nj
    C = reshape(reshape(A.value, n, n) * reshape(B.value, n, n), ni, nj, ni, nj)
    return ProductSpaceOperator(ProductSpace{S1, S2}, C)
end

# operator / field
function *(A::ProductSpaceOperator{ProductSpace{S1,S2}},
           u::Field{ProductSpace{S1, S2}})::Field{ProductSpace{S1, S2}} where {S1, S2}
    v = similar(u.value)
    for index in CartesianIndices(size(v))
        i, ii = index.I
        v[index] = sum(A.value[i, ii, k, kk]*u.value[k, kk] for k in range(S2), kk in range(S1))
    end
    return Field(ProductSpace{S1, S2}, v)
end

function *(u::Field{PS}, A::ProductSpaceOperator{PS})::ProductSpaceOperator{PS} where {PS}
    B = similar(A.value)
    for index in CartesianIndices(size(B))
        i, j, ii, jj = index.I
        B[index]     = u.value[i, j]*A.value[index]
    end
    return ProductSpaceOperator(PS, B)
end

function solve(A::ProductSpaceOperator{ProductSpace{S1, S2}}, u::Field{ProductSpace{S1, S2}})::Field{ProductSpace{S1, S2}} where {S1, S2}
    X =  vec(A) \ vec(u)
    return Field(ProductSpace{S1, S2}, shape(ProductSpace{S1, S2}, X))
end

# compute fields
function Field(PS::Type{ProductSpace{S1, S2}}, umap::Function)::Field{PS} where {S1, S2 <: Cardinal{Tag, N}} where {Tag, N}
    value = zeros(size(PS))
    for index in CartesianIndices(value)
        value[index] = umap(collocation(S1, index.I[1]),
                            collocation(S2, index.I[2]))
    end
    return Field(PS, value)
end

# compute fields at given U, V
function Field(PS::Type{ProductSpace{S1, S2}}, umap::Function,
               u::Field{ProductSpace{S1, S2}},
               v::Field{ProductSpace{S1, S2}})::Field{ProductSpace{S1, S2}} where {S1, S2 <: Cardinal{Tag,N}} where {Tag, N}
    value = zeros(spacetype(PS), size(PS))
    for index in range(PS)
        value[index] = umap(u.value[index], v.value[index])
    end
    return Field(PS, value)
end

# field / boundary
function +(u::Real, b::Boundary{ProductSpace{S1,S2}})::Field{ProductSpace{S1,S2}} where {S1,S2}
    return Field(ProductSpace{S1,S2}, u + b.value)
end

# scalar / boundary
function +(u::Field{ProductSpace{S1,S2}}, b::Boundary{ProductSpace{S1,S2}})::Field{ProductSpace{S1,S2}} where {S1,S2}
    return Field(ProductSpace{S1,S2}, u.value + b.value)
end

# Define multiplication rules to extract operators from J
function Base. *(A::ProductSpaceOperator{S}, u::Symbol)::ProductSpaceOperator{S} where {S}
    return A
end

function Base. *(v::Field{S}, u::Symbol)::ProductSpaceOperator{S} where {S<:ProductSpace{S1, S2}} where {S1, S2}
    return v*eye(S)
end

function Base. *(A::ProductSpaceOperator{S}, u::Int)::ProductSpaceOperator{S} where {S}
    @assert u == 0 "Right mutiplication with an Int is only allowed to represent \n action on a zero vector field"
    return (A-A)                 
end

function Base. *(v::Field{S}, u::Int)::ProductSpaceOperator{S} where {S<:ProductSpace{S1, S2}} where {S1, S2}
    @assert u == 0 "Right mutiplication with an Int is only allowed to represent \n action on a zero vector field"
    return (eye(S) - eye(S))                
end

function ⊙(A::ProductSpaceOperator{S}, B::ProductSpaceOperator{S})::ProductSpaceOperator{S} where {S}
    I = eye(S)
    return (I-B)*A + B*B
end

function ⊙(A::ProductSpaceOperator{S}, B::ProductSpaceOperator{S}, D::ProductSpaceOperator{S})::ProductSpaceOperator{S} where {S}
    I = eye(S)
    return (I-B)*A + B*D
end
