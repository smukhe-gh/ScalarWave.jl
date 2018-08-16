#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
#--------------------------------------------------------------------

import Base: *, +, -

abstract type Manifold end
abstract type Space <: Manifold end

struct Chebyshev{T} <: Space where T<:Real
    range::Range
    order::Int
end

struct Field{S<:Space, T<:Real} 
    space::S
    coefficents::Array{T,1}
end

struct Operator{Sin, Sout, T<:Real}
    space::Space
    value::Array{T,2}
end

function +{S}(u::Field{S}, v::Field{S})::Field{S}
    return Field(u.space, u.coefficents.+v.coefficents) 
end

function *(u::Field, v::Field)::Field
    # NOTE: Not sure if this will work in all cases
    @assert u.space === v.space
    return Field(u.space, u.coefficents.*v.coefficents) 
end

function derivative(S::Space)::Operator
    if typeof(S) == Chebyshev
        D = [chebd(i, j, S.order) for i in S.range, j in S.range]
        return Operator(S, D)
    else
        warn("Unknown space. Cannot compute the derivative operator")
        return Operator(S, eye(S.order+1))
end

function *(D::Operator, u::Field)::Field
    @assert D.space === u.space
    return Field(u.space, D.value*u.coefficents)
end

function *(D::Operator, W::Operator)::Operator
    @assert D.space === W.space
    return Field(u.space, D.value*W.value)
end

