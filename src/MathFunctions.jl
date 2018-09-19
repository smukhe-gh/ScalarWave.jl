#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Math functions on fields
#--------------------------------------------------------------------

import Base: sin, cos, sinpi, cospi,  exp, ^

sin(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, sin.(A.value))
cos(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, cos.(A.value))
sinpi(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, sinpi.(A.value))
cospi(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, cospi.(A.value))
exp(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Cardinal{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, exp.(A.value))

function ^(B::Field{ProductSpace{S1, S2}}, a::Int) where {S1, S2 <: Cardinal{Tag}} where {Tag}
    A = similar(B.value)
    for index in CartesianRange(size(A))
        A[index] = B.value[index].^a
    end
    return Field(ProductSpace{S1, S2}, A)
end
