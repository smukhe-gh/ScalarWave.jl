#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2019
# Test residuals away from axis
#--------------------------------------------------------------------

function computeSonUboundary(ϕ::Field{S})::Field{S} where {S<:Space{Tag}} where {Tag}
    v = Field(ϕ.space, x->x)
    DV = derivative(ϕ.space)
    I  = identity(ϕ.space)
    B  = incomingboundary(ϕ.space)
    b = 4*π*(DV*ϕ)^2
    L = 2*I + 2*v*DV 
    return solve(L ⊕ B, (I-B)*b)
end

function computeSonVboundary(ϕ::Field{S})::Field{S} where {S<:Space{Tag}} where {Tag}
    return -computeSonUboundary(ϕ)
end

function computeUboundary(ϕ::Field{S})::NTuple{3, Field{S}} where {S<:Space{Tag}} where {Tag}
    ψ = Field(ϕ.space, v->1)
    return (computeSonUboundary(ϕ), ψ, ϕ)
end

function computeVboundary(ϕ::Field{S})::NTuple{3, Field{S}} where {S<:Space{Tag}} where {Tag}
    ψ = Field(ϕ.space, u->1)
    return (computeSonVboundary(ϕ), ψ, ϕ)
end

PS = ProductSpace(ChebyshevGL{U, 3, Float64}(0, 1), 
                  ChebyshevGL{V, 3, Float64}(0, 1))

ϕ  = Field(PS, (u,v)->sin(u+v))
(s, ψ, ϕ) = combineUVboundary(computeUboundary(extractUboundary(ϕ, :incoming)), 
                              computeVboundary(extractVboundary(ϕ, :incoming)), :incoming)

display(s)
display(ψ)
display(ϕ)


