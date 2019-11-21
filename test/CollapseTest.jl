#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2019
# Simulate scalar field collpase on axis
#--------------------------------------------------------------------

using NLsolve, ForwardDiff
import LinearAlgebra.eigvals, LinearAlgebra.display, LinearAlgebra.cond

function C(a::Field{S}, η::Field{S}, ϕ::Field{S})::NTuple{2, Field{S}} where {S}
    C1 = DU*(DU*η) - (2/a)*(DU*a)*(DU*η) + (4pi*η)*(DU*ϕ)^2
    C2 = DV*(DV*η) - (2/a)*(DV*a)*(DV*η) + (4pi*η)*(DV*ϕ)^2
    return (C1, C2)
end

function diagnostics(a::Field{S}, η::Field{S}, ϕ::Field{S}) where {S}
    # Check if certain basic regularity conditions are being met on axis
    R1 = DU*η + DV*η
    R2 = DU*ϕ - DV*ϕ

    R3 = (1/2)*a - DV*η
    R4 = (1/2)*a + DU*η
    R5 = (1/4)*a^2 + (DU*η)*(DV*η)

    # println("--------------------------------------------------------------------------------")
    # display(A*(DV*η))
    # display(A*(DU*η))
    # display(A*((1/2)*a^2))
    # println("--------------------------------------------------------------------------------")
   
    @show L2(A*R1)
    @show L2(A*R2)
    @show L2(A*R3)
    @show L2(A*R4)
    @show L2(A*R5)

    # Check the domain of η. Would give us an idea of the spatial scales we're
    # trying to resolve
    # @show maximum(η)
    # @show minimum(η)

    # Check certain parity conditions on axis.  
    # Choose u = t - η, v = t + η. Reflect on this choice and what it means, 
    # but this let's us compute da/dη and dϕ/dη
end

function einstein(boundarydata::NTuple{3, Field{ProductSpace{S1, S2}}}, 
                  initialguess::NTuple{3, Field{ProductSpace{S1, S2}}})::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}

    function F(a::Field{S}, η::Field{S}, ϕ::Field{S})::NTuple{3, Field{S}} where {S}
        F1 = DU*(DV*a) - (1/a)*(DU*a)*(DV*a) + (a/η)*(DU*(DV*η)) + (4pi*a)*(DU*ϕ)*(DV*ϕ)
        F2 = DU*(DV*η) + (1/η)*(DU*η)*(DV*η) + (1/4)*(1/η)*(a^2)
        F3 = DU*(DV*ϕ) + (1/η)*(DU*η)*(DV*ϕ) + (1/η)*(DV*η)*(DU*ϕ)

        return (F1, F2, F3)

        F1onAxis = DU*DV*a - (1/a)*(DU*a)*(DV*a) + (4pi*a)*(DU*ϕ)*(DV*ϕ)
        F1onAxis = (DU*η)*(DV*η) + (1/4)*a^2 
        # F2onAxis = DU*DV*η 
        # F3onAxis = DU*DV*ϕ 
        F3onAxis = DU*ϕ - DV*ϕ 

        if typeof(a.value) <: Array{Float64, 2}
            diagnostics(a, η, ϕ)
        end

        return (mix!(mix!(F1, A, F1onAxis), B, a-bnda), 
                mix!(mix!(F2, A, F2onAxis), B, η-bndη),
                mix!(mix!(F3, A, F3onAxis), B, ϕ-bndϕ))
    end

    function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
        fvec[:] = reshapeFromTuple(F(reshapeToTuple(PS, x)...))
    end

    
    function j!(jvec::Array{T,2}, x::Array{T,1}) where {T}
        jvec[:, :] = ForwardDiff.jacobian(f!, similar(x), x)
        display(jvec)
        exit()
    end

    (a0, η0, ϕ0) = initialguess
    (bnda, bndη, bndϕ) = combineUVboundary.(computeUboundary((a0, η0, ϕ0)), 
                                            computeVboundary((a0, η0, ϕ0)), :incoming)

    @show lineconstraint(extractUboundary.((bnda, bndη, bndϕ), :incoming)...)
    @show lineconstraint(extractVboundary.((bnda, bndη, bndϕ), :incoming)...)
    @show L2.(C(a0, η0, ϕ0))
    @show L2.(F(a0, η0, ϕ0))

    # solved = reshapeToTuple(PS, nlsolve(f!, j!,  reshapeFromTuple((a0, η0, ϕ0)); method=:trust_region, show_trace=true, ftol=1e-10, iterations=120).zero)

    solved = reshapeToTuple(PS, nlsolve(f!, reshapeFromTuple((a0, η0, ϕ0)); method=:trust_region, autodiff=:forward,
                                            show_trace=true, ftol=1e-10, iterations=120).zero)

    @show L2.(C(solved...))

    return solved
end

PS = ProductSpace(ChebyshevGL{U, 12, Float64}(0, 1), 
                  ChebyshevGL{V, 12, Float64}(1.1, 2.1))

DU, DV = derivative(PS)
B = incomingboundary(PS)
A = axisboundary(PS)
I = identity(PS)

p  = 1e-1
a0 = Field(PS, (u,v)->1)
η0 = Field(PS, (u,v)->(v-u)/2)
ϕ0 = Field(PS, (u,v)->p*exp(-(v-2)^2)
                    + p*exp(-(u-2)^2))

(asol, ηsol, ϕsol) = einstein((B*a0, B*η0, B*ϕ0), (a0, η0, ϕ0))

