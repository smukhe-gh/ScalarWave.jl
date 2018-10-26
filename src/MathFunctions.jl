#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Math functions on constructed datatypes 
#--------------------------------------------------------------------

import Base: sqrt, abs, sin, cos, sinpi, cospi,  exp, ^, log
import Base: maximum, minimum

sqrt(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Space{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, sqrt.(A.value))
abs(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Space{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, abs.(A.value))
sin(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Space{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, sin.(A.value))
cos(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Space{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, cos.(A.value))
sinpi(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Space{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, sinpi.(A.value))
cospi(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Space{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, cospi.(A.value))
exp(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Space{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, exp.(A.value))
log(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Space{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, log.(A.value))
^(B::Field{ProductSpace{S1, S2}}, a::Int) where {S1, S2 <: Space{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, B.value.^a)

maximum(u::Field{ProductSpace{S1, S2}}) where {S1, S2} = maximum(u.value)
minimum(u::Field{ProductSpace{S1, S2}}) where {S1, S2}  = minimum(u.value)

import LinearAlgebra: eigvals, cond

function eigvals(A::ProductSpaceOperator{ProductSpace{S1, S2}}) where {S1, S2}
    return eigvals(vec(A))
end

function cond(A::ProductSpaceOperator{ProductSpace{S1, S2}}) where {S1, S2}
    return cond(vec(A))
end

