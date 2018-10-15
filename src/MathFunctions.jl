#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Math functions on fields
#--------------------------------------------------------------------

import Base: sin, cos, sinpi, cospi,  exp, ^, log, sum

sin(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, sin.(A.value))
cos(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, cos.(A.value))
sinpi(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, sinpi.(A.value))
cospi(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, cospi.(A.value))
exp(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, exp.(A.value))
log(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, log.(A.value))

function ^(B::Field{ProductSpace{S1, S2}}, a::Int) where {S1, S2 <: Cardinal{Tag}} where {Tag}
    A = similar(B.value)
    for index in CartesianRange(size(A))
        A[index] = B.value[index].^a
    end
    return Field(ProductSpace{S1, S2}, A)
end

# Adding this here to to a lack of a better place
# Compute L2 error norm

sum(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, sum(A.value))

function L2norm(A::Field{ProductSpace{S1, S2}})::Float64 where {S1, S2 <: Cardinal{Tag}} where {Tag}
    W = [chebw(i, order(S2))*chebw(j,order(S1)) for i in 1:len(S2), j in 1:len(S1)]
    𝕎 = ProductSpaceOperator(ProductSpace{S1, S2}, diagm(vec(W)))
    return sqrt(sum(𝕎 *(A^2)).value)
end

