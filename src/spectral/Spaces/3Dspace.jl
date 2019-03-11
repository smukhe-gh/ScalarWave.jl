#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Define basic math and array operations for 3D spaces
#--------------------------------------------------------------------

import Base: +, -, *, /, size, range, reshape, zero, similar, maximum, minimum

# dimensions and shape
order(PS::Type{ProductSpace{S1, S2, S3}}) where {S1, S2, S3} = (order(S3), order(S2), order(S1))
dim(PS::Type{ProductSpace{S1, S2, S3}}) where {S1, S2, S3} = dim(S1) + dim(S2) + dim(S3)
range(PS::Type{ProductSpace{S1, S2, S3}}) where {S1, S2, S3} = CartesianIndices((length(S3), length(S2), length(S1)))
size(PS::Type{ProductSpace{S1, S2, S3}}) where {S1, S2, S3} = (length(S3), length(S2), length(S1))

# typeof 
spacetype(::Type{ProductSpace{S1, S2, S3}}) where {S1, S2, S3 <: GaussLobatto{Tag,N}} where {Tag, N} = Float64
spacetype(::Type{ProductSpace{S1, S2, S3}}) where {S1, S2, S3 <: Taylor{Tag,N}} where {Tag, N} = Rational{BigInt}

# maximum and minimum
maximum(::Type{ProductSpace{S1, S2, S3}}) where {S1, S2, S3 <: GaussLobatto{Tag,N}} where {Tag, N} = (maximum(S1), maximum(S2), maximum(S3))
minimum(::Type{ProductSpace{S1, S2, S3}}) where {S1, S2, S3 <: GaussLobatto{Tag,N}} where {Tag, N} = (minimum(S1), minimum(S2), maximum(S3))

# identity element
zero(::Type{ProductSpace{S1, S2, S3}}) where {S1, S2, S3} = Field(ProductSpace{S1, S2, S3}, (u,v)->0)
zero(::Type{Null}, ::Type{ProductSpace{S1, S2, S3}}) where {S1, S2, S3} = 0*boundary(Null, ProductSpace{S1, S2, S3})

# unary operators
-(A::Field{ProductSpace{S1, S2, S3}}) where {S1, S2, S3 <: Cardinal{Tag, N}} where {Tag,N} = Field(ProductSpace{S1, S2, S3}, -A.value)
+(A::Field{ProductSpace{S1, S2, S3}}) where {S1, S2, S3 <: Cardinal{Tag, N}} where {Tag,N} = Field(ProductSpace{S1, S2, S3}, +A.value)

# field / field
+(A::Field{ProductSpace{S1, S2, S3}},
  B::Field{ProductSpace{S1, S2, S3}}) where {S1, S2, S3 <: Cardinal{Tag, N}} where {Tag,N} = Field(ProductSpace{S1, S2, S3}, A.value .+ B.value)
-(A::Field{ProductSpace{S1, S2, S3}},
  B::Field{ProductSpace{S1, S2, S3}}) where {S1, S2, S3 <: Cardinal{Tag, N}} where {Tag,N} = Field(ProductSpace{S1, S2, S3}, A.value .- B.value)
*(A::Field{ProductSpace{S1, S2, S3}},
  B::Field{ProductSpace{S1, S2, S3}}) where {S1, S2, S3 <: Cardinal{Tag, N}} where {Tag,N} = Field(ProductSpace{S1, S2, S3}, A.value .* B.value)
/(A::Field{ProductSpace{S1, S2, S3}},
  B::Field{ProductSpace{S1, S2, S3}}) where {S1, S2, S3 <: Cardinal{Tag, N}} where {Tag,N} = Field(ProductSpace{S1, S2, S3}, A.value ./ B.value)

# scalar / field
+(a::T,
  B::Field{ProductSpace{S1, S2, S3}}) where {S1, S2, S3 <: Cardinal{Tag, N}, T<:Real} where {Tag,N} = Field(ProductSpace{S1, S2, S3}, a .+ B.value)
-(a::T,
  B::Field{ProductSpace{S1, S2, S3}}) where {S1, S2, S3 <: Cardinal{Tag, N}, T<:Real} where {Tag,N} = Field(ProductSpace{S1, S2, S3}, a .- B.value)
*(a::T,
  B::Field{ProductSpace{S1, S2, S3}}) where {S1, S2, S3 <: Cardinal{Tag, N}, T<:Real} where {Tag,N} = Field(ProductSpace{S1, S2, S3}, a .* B.value)
/(a::T,
  B::Field{ProductSpace{S1, S2, S3}}) where {S1, S2, S3 <: Cardinal{Tag, N}, T<:Real} where {Tag,N} = Field(ProductSpace{S1, S2, S3}, a ./ B.value)
/(B::Field{ProductSpace{S1, S2, S3}}, a::T) where {S1, S2, S3 <: Cardinal{Tag, N}, T<:Real} where {Tag,N} = Field(ProductSpace{S1, S2, S3}, B.value ./ a)

+(a::T,
  B::Field{ProductSpace{S1, S2, S3}}) where {S1, S2, S3 <: Space{Tag}, T<:Complex} where {Tag} = Field(ProductSpace{S1, S2, S3}, a .+ B.value)
-(a::T,
  B::Field{ProductSpace{S1, S2, S3}}) where {S1, S2, S3 <: Space{Tag}, T<:Complex} where {Tag} = Field(ProductSpace{S1, S2, S3}, a .- B.value)
*(a::T,
  B::Field{ProductSpace{S1, S2, S3}}) where {S1, S2, S3 <: Space{Tag}, T<:Complex} where {Tag} = Field(ProductSpace{S1, S2, S3}, a .* B.value)

function *(A::ProductSpaceOperator{ProductSpace{S1,S2,S3}},
           B::ProductSpaceOperator{ProductSpace{S1,S2,S3}})::ProductSpaceOperator{ProductSpace{S1, S2, S3}} where {S1, S2, S3}
    nh = length(S3)
    ni = length(S2)
    nj = length(S1)
    n = nh * ni * nj
    C = reshape(reshape(A.value, n, n) * reshape(B.value, n, n), nh, ni, nj, nh, ni, nj)
    return ProductSpaceOperator(ProductSpace{S1, S2, S3}, C)
end

# operator / field
function *(A::ProductSpaceOperator{ProductSpace{S1,S2,S3}},
           u::Field{ProductSpace{S1, S2, S3}})::Field{ProductSpace{S1, S2, S3}} where {S1, S2, S3}
    v = similar(u.value)
    for index in CartesianIndices(size(v))
        i, j, k = index.I
        v[index] = sum(A.value[i, j, k, l, m, n]*u.value[l, m, n] for l in range(S3), m in range(S2), n in range(S1))
    end
    return Field(ProductSpace{S1, S2, S3}, v)
end

function *(u::Field{ProductSpace{S1, S2, S3}}, A::ProductSpaceOperator{ProductSpace{S1, S2, S3}})::ProductSpaceOperator{ProductSpace{S1, S2, S3}} where {S1, S2, S3}
    B = similar(A.value)
    for index in CartesianIndices(size(B))
        i, j, k, l, m, n = index.I
        B[index]     = u.value[i, j, k]*A.value[index]
    end
    return ProductSpaceOperator(ProductSpace{S1, S2, S3}, B)
end

# TODO: Check if solves correctly.
function solve(A::ProductSpaceOperator{ProductSpace{S1, S2, S3}}, u::Field{ProductSpace{S1, S2, S3}})::Field{ProductSpace{S1, S2, S3}} where {S1, S2, S3}
    X =  reshape(A) \ reshape(u)
    return reshape(ProductSpace{S1, S2, S3}, X)
end

# compute fields
function Field(PS::Type{ProductSpace{S1, S2, S3}}, umap::Function)::Field{PS} where {S1, S2, S3 <: Cardinal{Tag, N}} where {Tag, N}
    value = zeros(spacetype(PS), size(PS))
    for index in range(PS)
        value[index] = umap(collocation(S3, index.I[1]),
                            collocation(S2, index.I[2]),
                            collocation(S1, index.I[3]))
    end
    return Field(PS, value)
end

# Kronecker product of three 1D operators 
function ⦼(A::Operator{S1}, B::Operator{S2}, C::Operator{S3})::ProductSpaceOperator{ProductSpace{S1, S2, S3}} where {S1, S2, S3<: Cardinal{Tag, N}} where {Tag, N}
    ABC = zeros(spacetype(A.space), length(S3), length(S2), length(S1), 
                                    length(S3), length(S2), length(S1))
    for index in CartesianIndices(size(ABC))
        i ,j , k, ii ,jj, kk = index.I
        ABC[index]    = A.value[j,jj]*B.value[i,ii]*C.value[k,kk]
    end
    return ProductSpaceOperator(ProductSpace{S1, S2, S3}, ABC)
end

# derivative operators
function derivative(PS::Type{ProductSpace{S1, S2, S3}}) where {S1, S2, S3}
    return (reshape(ProductSpace{S1, S2, S3}, kron(derivative(S1).value, eye(S2).value, eye(S3).value)),
            reshape(ProductSpace{S1, S2, S3}, kron(eye(S1).value, derivative(S2).value, eye(S3).value)),
            reshape(ProductSpace{S1, S2, S3}, kron(eye(S1).value, eye(S2).value, derivative(S3).value)))
end

# create the identity operator
function eye(PS::Type{ProductSpace{S1, S2, S3}}) where {S1, S2, S3}
    return reshape(ProductSpace{S1, S2, S3}, kron(eye(S1).value,  eye(S2).value, eye(S3).value))
end

# Null and Spatial boundary operators
function boundary(::Type{Spacelike}, PS::Type{ProductSpace{S1, S2, S3}})::ProductSpaceOperator{ProductSpace{S1, S2, S3}} where {S1, S2, S3 <: Cardinal{Tag,N}}  where {Tag, N}
    B = zeros(spacetype(PS), length(S3), length(S2), length(S1))
    B[1, :, :] = B[end, :, :] .= convert(spacetype(PS), 1)
    B[:, :, 1] = B[:, :, end] .= convert(spacetype(PS), 1)
    B[:, 1, :] = B[:, end, :] .= convert(spacetype(PS), 1)
    return ProductSpaceOperator(ProductSpace{S1, S2, S3}, reshape(diagm(0=>vec(B)), (length(S3), length(S2), length(S1), 
                                                                                     length(S3), length(S2), length(S1))))
end

function boundary(::Type{Null}, PS::Type{ProductSpace{S1, S2, S3}})::ProductSpaceOperator{ProductSpace{S1, S2, S3}} where {S1, S2, S3 <: Cardinal{Tag,N}}  where {Tag, N}
    B = zeros(spacetype(PS), length(S3), length(S2), length(S1))
    B[1, :, :] .= convert(spacetype(PS), 1)
    B[:, :, 1] .= convert(spacetype(PS), 1)
    B[:, 1, :] .= convert(spacetype(PS), 1)
    return ProductSpaceOperator(ProductSpace{S1, S2, S3}, reshape(diagm(0=>vec(B)), (length(S3), length(S2), length(S1), 
                                                                                     length(S3), length(S2), length(S1))))
end

# shape and reshape operators
function reshape(u::Field{ProductSpace{S1, S2, S3}})::Array{eltype(u.value),1} where {S1, S2, S3}
    return reshape(u.value, prod(size(ProductSpace{S1, S2, S3})))
end

function reshape(S::Type{ProductSpace{S1, S2, S3}}, A::Array{Float64,1})::Field{ProductSpace{S1, S2, S3}} where {S1, S2, S3}
    return Field(ProductSpace{S1, S2, S3}, reshape(A, (length(S3), length(S2), length(S1))))
end

function reshape(A::ProductSpaceOperator{ProductSpace{S1, S2, S3}})::Array{eltype(A.value),2} where {S1, S2, S3}
    return reshape(A.value, (prod(size(ProductSpace{S1, S2, S3})), prod(size(ProductSpace{S1, S2, S3}))))
end

function reshape(S::Type{ProductSpace{S1, S2, S3}}, A::Array{Float64,2})::ProductSpaceOperator{ProductSpace{S1, S2, S3}} where {S1, S2, S3}
    return ProductSpaceOperator(ProductSpace{S1, S2, S3}, reshape(A, (length(S3), length(S2), length(S1),
                                                                      length(S3), length(S2), length(S1))))
end
