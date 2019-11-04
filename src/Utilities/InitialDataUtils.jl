#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 09-2019
# Utilities for initial data solver
#--------------------------------------------------------------------

export solveη, lineconstraint, computeUboundary, computeVboundary

function solveη(a::Field{S}, η::Field{S}, ϕ::Field{S})::Field{S} where {S<:Space{Tag}} where {Tag}
    D = derivative(a.space)
    I = identity(a.space)
    B = incomingboundary(a.space) + outgoingboundary(a.space)
    A = D*D - (2/a)*(D*a)*D + 4pi*(D*ϕ)^2*I
    return solve(A ⊕ B, B*η)
end

function lineconstraint(a::Field{S}, η::Field{S}, ϕ::Field{S})::Number where {S<:Space{Tag}} where {Tag}
    D = derivative(a.space)
    h = D*(D*η) - (2/a)*(D*a)*(D*η) + (4pi*η)*(D*ϕ)^2
    return L2(h)
end

function computeUboundary(u::NTuple{3, Field{ProductSpace{S1, S2}}})::NTuple{3, Field{S2}} where {S1, S2}
    return (extractUboundary(u[1], :incoming), solveη(extractUboundary.(u, :incoming)...), extractUboundary(u[3], :incoming))
end

function computeVboundary(u::NTuple{3, Field{ProductSpace{S1, S2}}})::NTuple{3, Field{S1}} where {S1, S2}
    return (extractVboundary(u[1], :incoming), solveη(extractVboundary.(u, :incoming)...), extractVboundary(u[3], :incoming))
end

