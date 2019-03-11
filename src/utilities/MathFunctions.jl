#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Math functions on constructed datatypes 
#--------------------------------------------------------------------

import Base: sqrt, abs, sin, cos, sinpi, cospi,  exp, ^, log, log10
import Base: maximum, minimum

sqrt(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Space{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, sqrt.(A.value))
abs(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Space{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, abs.(A.value))
sin(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Space{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, sin.(A.value))
cos(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Space{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, cos.(A.value))
sinpi(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Space{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, sinpi.(A.value))
cospi(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Space{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, cospi.(A.value))
exp(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Space{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, exp.(A.value))
log(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Space{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, log.(A.value))
log10(A::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Space{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, log10.(A.value))
^(B::Field{ProductSpace{S1, S2}}, a::Int) where {S1, S2 <: Space{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, (B.value).^a)
^(B::Field{S}, a::Int) where {S <: Cardinal{Tag}} where {Tag} = Field(S, (B.value).^a)
inv(B::Field{ProductSpace{S1, S2}}) where {S1, S2 <: Space{Tag}} where {Tag} = Field(ProductSpace{S1, S2}, 1 ./(B.value))

maximum(u::Field{ProductSpace{S1, S2}}) where {S1, S2} = maximum(u.value)
minimum(u::Field{ProductSpace{S1, S2}}) where {S1, S2}  = minimum(u.value)

import LinearAlgebra: eigvals, cond

function eigvals(A::ProductSpaceOperator{ProductSpace{S1, S2}}) where {S1, S2}
    return eigvals(vec(A))
end

function cond(A::ProductSpaceOperator{ProductSpace{S1, S2}}) where {S1, S2}
    return cond(vec(A))
end

function cond(A::ProductSpaceOperator{ProductSpace{S1, S2, S3}}) where {S1, S2, S3}
    return cond(reshape(A))
end

abs(A::Field{ProductSpace{S1, S2, S3}}) where {S1, S2, S3<: Space{Tag}} where {Tag} = Field(ProductSpace{S1, S2, S3}, abs.(A.value))
maximum(u::Field{ProductSpace{S1, S2, S3}}) where {S1, S2, S3} = maximum(u.value)

# FIXME: Remove repetition for fields in 1D space
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


