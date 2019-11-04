#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2019
# Simulate Minkowski spacetime on axis
#--------------------------------------------------------------------

using NLsolve, ForwardDiff
import LinearAlgebra.eigvals, LinearAlgebra.display, LinearAlgebra.cond

function C(a::Field{S}, η::Field{S}, ϕ::Field{S})::NTuple{2, Field{S}} where {S}
    C1 = DU*(DU*η) - (2/a)*(DU*a)*(DU*η) + (4pi*η)*(DU*ϕ)^2
    C2 = DV*(DV*η) - (2/a)*(DV*a)*(DV*η) + (4pi*η)*(DV*ϕ)^2
    return (C1, C2)
end

function einstein(boundarydata::NTuple{3, Field{ProductSpace{S1, S2}}}, 
                  initialguess::NTuple{3, Field{ProductSpace{S1, S2}}})::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}

    function F(a::Field{S}, η::Field{S}, ϕ::Field{S})::NTuple{3, Field{S}} where {S}
        F1 = DU*(DV*a) - (1/a)*(DU*a)*(DV*a) + (a/η)*(DU*(DV*η)) + (4pi*a)*(DU*ϕ)*(DV*ϕ)
        F2 = DU*(DV*η) + (1/η)*(DU*η)*(DV*η) + (1/4)*(1/η)*(a^2)
        F3 = DU*(DV*ϕ) + (1/η)*(DU*η)*(DV*ϕ) + (1/η)*(DV*η)*(DU*ϕ)

        F1onAxis = DU*DV*a - (1/a)*(DU*a)*(DV*a) + (4pi*a)*(DU*ϕ)*(DV*ϕ)
        F2onAxis = DU*DV*η 
        F3onAxis = DU*DV*ϕ 

        return (mix!(F1, B, a-bnda), 
                mix!(F2, B, η-bndη),
                mix!(F3, B, ϕ-bndϕ))

        # return (mix!(mix!(F1, A, F1onAxis), B, a-bnda), 
                # mix!(mix!(F2, A, F2onAxis), B, η-bndη),
                # mix!(mix!(F3, A, F3onAxis), B, ϕ-bndϕ))
    end

    function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
        fvec[:] = reshapeFromTuple(F(reshapeToTuple(PS, x)...))
    end

    
    function j!(jvec::Array{T,2}, x::Array{T,1}) where {T}
        jvec[:, :] = ForwardDiff.jacobian(f!, similar(x), x)
    end

    (a0, η0, ϕ0) = initialguess
    (bnda, bndη, bndϕ) = combineUVboundary.(computeUboundary((a0, η0, ϕ0)), 
                                            computeVboundary((a0, η0, ϕ0)), :incoming)

    @show lineconstraint(extractUboundary.((bnda, bndη, bndϕ), :incoming)...)
    @show lineconstraint(extractVboundary.((bnda, bndη, bndϕ), :incoming)...)
    @show L2.(C(a0, η0, ϕ0))
    @show L2.(F(a0, η0, ϕ0))

    solved = reshapeToTuple(PS, nlsolve(f!, reshapeFromTuple((a0, η0, ϕ0)); method=:trust_region, autodiff=:forward,
                                            show_trace=true, ftol=1e-10, iterations=100).zero)

    @show L2.(C(solved...))
    @show L2.(F(solved...))

    return solved
end

PS = ProductSpace(ChebyshevGL{U, 30, Float64}(0, 1), 
                  ChebyshevGL{V, 30, Float64}(1.1, 2.1))

DU, DV = derivative(PS)
B = incomingboundary(PS)
A = axisboundary(PS)
I = identity(PS)

p  = 1e-2
a0 = Field(PS, (u,v)->1)
η0 = Field(PS, (u,v)->(v-u)/2)
ϕ0 = Field(PS, (u,v)->p*exp(-(v-1)^2))

(asol, ηsol, ϕsol) = einstein((B*a0, B*η0, B*ϕ0), (a0, η0, ϕ0))

