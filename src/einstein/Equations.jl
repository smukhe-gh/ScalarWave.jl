#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# Bad choice for a file name. This holds the expressions for
# linear and non-operators
# There are 4 equations: two hyperbolic in nature (H1, H2)
# and two elliptic in nature (E3, E4). 
# See Gundlach and Pullin 1997
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Compute residuals for non-linear equations
#--------------------------------------------------------------------

function H1(f::Field{S}, r::Field{S}, ϕ::Field{S})::Field{S} where {S<:Cardinal{Tag}} where {Tag}
    DV, DU = derivative(S)
    return DV*DU*r + ((DU*r)*(DV*r))/r - f/r
end

function H2(f::Field{S}, r::Field{S}, ϕ::Field{S})::Field{S} where {S<:Cardinal{Tag}} where {Tag}
    DV, DU = derivative(S)
    return DV*DU*log(f) + (2/r)*DV*DU*r + 2*(DU*ϕ)*(DV*ϕ)
end

function E1(f::Field{S}, r::Field{S}, ϕ::Field{S})::Field{S} where {S<:Cardinal{Tag}} where {Tag}
    DV, DU = derivative(S)
    return DU*DU*r - ((DU*f)/f)*DU*r + (DU*ϕ)^2
end

function E2(f::Field{S}, r::Field{S}, ϕ::Field{S})::Field{S} where {S<:Cardinal{Tag}} where {Tag}
    DV, DU = derivative(S)
    return DV*DV*r - ((DV*f)/f)*DV*r + (DV*ϕ)^2
end

#--------------------------------------------------------------------
# Compute linear operators 
# TODO: Compute identity operator
#--------------------------------------------------------------------

function linearH1(f::Field{S}, r::Field{S}, ϕ::Field{S})::ProductSpaceOperator{S} where {S<:Cardinal{Tag}} where {Tag}
end

function linearH2(f::Field{S}, r::Field{S}, ϕ::Field{S})::ProductSpaceOperator{S} where {S<:Cardinal{Tag}} where {Tag}
end

function linearE1(f::Field{S}, r::Field{S}, ϕ::Field{S})::Field{S} where {S<:Cardinal{Tag}} where {Tag}
end

function linearE2(f::Field{S}, r::Field{S}, ϕ::Field{S})::Field{S} where {S<:Cardinal{Tag}} where {Tag}
end
