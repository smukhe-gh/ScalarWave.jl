#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Math functions on constructed datatypes 
#--------------------------------------------------------------------

import Base: sqrt, abs, sin, cos, sinpi, cospi,  exp, ^, log
import Base: maximum, minimum

sqrt(A::Field{S}) where {S <: Space{Tag}} where {Tag} = Field(S, sqrt.(A.value))
abs(A::Field{S}) where {S <: Space{Tag}} where {Tag} = Field(S, abs.(A.value))
sin(A::Field{S}) where {S <: Space{Tag}} where {Tag} = Field(S, sin.(A.value))
cos(A::Field{S}) where {S <: Space{Tag}} where {Tag} = Field(S, cos.(A.value))
sinpi(A::Field{S}) where {S <: Space{Tag}} where {Tag} = Field(S, sinpi.(A.value))
cospi(A::Field{S}) where {S <: Space{Tag}} where {Tag} = Field(S, cospi.(A.value))
exp(A::Field{S}) where {S <: Space{Tag}} where {Tag} = Field(S, exp.(A.value))
log(A::Field{S}) where {S <: Space{Tag}} where {Tag} = Field(S, log.(A.value))
^(B::Field{S}, a::Int) where {S <: Space{Tag}} where {Tag} = Field(S, (B.value).^a)
inv(B::Field{S}) where {S <: Space{Tag}} where {Tag} = Field(S, 1 ./(B.value))

maximum(u::Field{S}) where {S <: Space{Tag}} where {Tag} = maximum(u.value)
minimum(u::Field{S}) where {S <: Space{Tag}} where {Tag} = minimum(u.value)

import LinearAlgebra: eigvals, cond

function eigvals(A::ProductSpaceOperator{ProductSpace{S1, S2}}) where {S1, S2}
    return eigvals(vec(A))
end

function cond(A::ProductSpaceOperator{ProductSpace{S1, S2}}) where {S1, S2}
    return cond(vec(A))
end

