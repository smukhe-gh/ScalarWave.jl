#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2019
# Simulate scalar field collpase on axis
#--------------------------------------------------------------------

using NLsolve, ForwardDiff
import LinearAlgebra.eigvals, LinearAlgebra.display, LinearAlgebra.cond

function operators(initialguess::NTuple{3, Field{ProductSpace{S1, S2}}})::NTuple{4, Operator{ProductSpace{S1, S2}}} where {S1, S2}
    PS = initialguess[1].space
    DU, DV = derivative(PS)
    A = axisboundary(PS)
    B = incomingboundary(PS)
    return (DU, DV, A, B)
end

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

function main(boundarydata::NTuple{3, Field{ProductSpace{S1, S2}}}, 
              initialguess::NTuple{3, Field{ProductSpace{S1, S2}}})::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}

    function C(s::Field{S}, ψ::Field{S}, ϕ::Field{S})::NTuple{2, Field{S}} where {S}
        C1 = (4*π*(DU*ϕ)^2 - 4*(r^2)*ψ*(DU*s)*(DU*ψ) 
            - 2*s*ψ*(DU*r)*(ψ*(DU*r) + 2*r*(DU*ψ)) 
            + (ψ^2)*(DU*(DU*r)) - 2*r*(ψ^2*(DU*r)*(DU*s) + 3*(DU*ψ)^2 - ψ*(DU*(DU*ψ))))
        C2 = (4*π*(DV*ϕ)^2 - 4*(r^2)*ψ*(DV*s)*(DV*ψ) 
            - 2*s*ψ*(DV*r)*(ψ*(DV*r) + 2*r*(DV*ψ)) 
            + (ψ^2)*(DV*(DV*r)) - 2*r*(ψ^2*(DV*r)*(DV*s) + 3*(DV*ψ)^2 - ψ*(DV*(DV*ψ))))
        return (C1, C2)
    end

    function F(s::Field{S}, ψ::Field{S}, ϕ::Field{S})::NTuple{3, Field{S}} where {S}
        F1 = ((1/r)*(ψ^2)*exp(2*r*s) + (1/r)*(ψ^2)*(DU*r)*(DV*r) 
                + 4*ψ*(DV*r)*(DU*ψ) + 4*ψ*(DU*r)*(DV*ψ) 
                + 6*r*(DU*ψ)*(DV*ψ) + (ψ^2)*(DU*(DV*r)) + 2*r*ψ*(DU*(DV*ψ)))
        F2 = (2*(DV*r)*(DU*s) + 2*(DU*r)*(DV*s)
                + 4*π*(DU*ϕ)*(DV*ϕ) + (4/(r*ψ))*((DV*r)*(DU*ψ) + (DU*r)*(DV*ψ))
                + (2/r)*(DU*(DV*r)) + 2*r*(DU*(DV*s)) + (8/ψ)*(DU*(DV*ψ)))
        F3 = ((1/r)*(DV*r)*(DU*ϕ) + (1/r)*(DU*r)*(DV*ϕ) 
                    + (2/ψ)*(DV*ϕ)*(DU*ψ) + (2/ψ)*(DU*ϕ)*(DV*ψ) + (DU*(DV*ϕ)))

        F1onAxis = s
        F2onAxis = DU*ψ - DV*ψ 
        F3onAxis = DU*ϕ - DV*ϕ 

        return (mix!(mix!(F1, A, F1onAxis), B, s-bnds), 
                mix!(mix!(F2, A, F2onAxis), B, ψ-bndψ),
                mix!(mix!(F3, A, F3onAxis), B, ϕ-bndϕ))
    end

    function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
        fvec[:] = reshapeFromTuple(F(reshapeToTuple(PS, x)...))
    end
    
    (DU, DV, A, B) = operators(initialguess)
    (s0, ψ0, ϕ0) = initialguess
    (bnds, bndψ, bndϕ) = boundarydata 

    @show L2.(C(initialguess...))
    solved = reshapeToTuple(PS, nlsolve(f!, reshapeFromTuple((s0, ψ0, ϕ0)); method=:trust_region, autodiff=:forward,
                                            show_trace=true, ftol=1e-10, iterations=120).zero)
    @show L2.(C(solved...))
    return solved
end

PS = ProductSpace(ChebyshevGL{U, 11, Float64}(0, 1), 
                  ChebyshevGL{V, 11, Float64}(0, 1))

ϕ  = Field(PS, (u,v)->sin(u+v))
boundarydata = combineUVboundary(computeUboundary(extractUboundary(ϕ, :incoming)), 
                                 computeVboundary(extractVboundary(ϕ, :incoming)), :incoming)

r = Field(PS, (u,v)->(v-u))
s = r
ψ = Field(PS, (u,v)->1)
initialguess = (s, ψ, ϕ)

main(boundarydata, initialguess)

