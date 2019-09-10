#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 07-2019
# Test functions inside functions
#--------------------------------------------------------------------

using NLsolve, PyPlot

function solveR(a::Field{S}, r::Field{S}, ϕ::Field{S})::Field{S} where {S<:Space{Tag}} where {Tag}
    D = derivative(ϕ.space)
    I = identity(ϕ.space)
    B = incomingboundary(ϕ.space) + outgoingboundary(ϕ.space)
    A = D*D - (2/a)*(D*a)*D + 4pi*(D*ϕ)^2*I
    return solve(A ⊕ B, B*r)
end

function initialdatasolver(a::Field{S}, r::Field{S}, ϕ::Field{S})::NTuple{3, Field{S}} where {S}
    rU = solveR(extractUboundary(a), extractUboundary(r), extractUboundary(ϕ))
    rV = solveR(extractVboundary(a), extractVboundary(r), extractVboundary(ϕ))
    rs  = combineUVboundary(rU, rV)
    B  = incomingboundary(a.space)
    return (B*a, B*rs, B*ϕ)
end

function nonlinearsolver(PS::ProductSpace{S1, S2})::NTuple{5, Field{ProductSpace{S1, S2}}} where {S1, S2 <: ChebyshevGL{Tag, N, T}} where {Tag, N, T}

    function C(a::Field{S}, r::Field{S}, ϕ::Field{S})::NTuple{2, Field{S}} where {S}
        C1 = DU*(DU*r) - (2/a)*(DU*a)*(DU*r) + (4pi*r)*(DU*ϕ)^2
        C2 = DV*(DV*r) - (2/a)*(DV*a)*(DV*r) + (4pi*r)*(DV*ϕ)^2
        return (C1, C2)
    end
    
    function F(a::Field{S}, r::Field{S}, ϕ::Field{S})::NTuple{3, Field{S}} where {S}
        F1 = r*(DU*(DV*ϕ)) + (DU*r)*(DV*ϕ) + (DV*r)*(DU*ϕ)
        F2 = r*(DU*(DV*r)) + (DU*r)*(DV*r) + (1/4)*(a^2)
        F3 = (1/a)*(DU*(DV*a)) - (1/a^2)*(DU*a)*(DV*a) + (1/r)*(DU*(DV*r)) + 4pi*(DU*ϕ)*(DV*ϕ)
        return (F1, F2, F3)
    end
    
    function residual(a::Field{S}, r::Field{S}, ϕ::Field{S})::NTuple{3, Field{S}} where {S}

        function BCs(F1::Field{S}, F2::Field{S}, F3::Field{S})::NTuple{3, Field{S}} where {S}
            BF1 = (I-B)*F1 + B*(a-bnda) 
            BF2 = (I-B)*F2 + B*(r-bndr) 
            BF3 = (I-B)*F3 + B*(ϕ-bndϕ) 
            return (BF1, BF2, BF3)
        end

        return BCs(F(a, r, ϕ)...)
    end

    function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
        fvec[:] = reshapeFromTuple(residual(reshapeToTuple(PS, x)...))
    end

    # Compute common operators
    DU, DV = derivative(PS)
    B = incomingboundary(PS)
    I = identity(PS)

    # Schwarzschild spacetime 
    # M  = T(1.0)
    # r0 = Field(PS, (u,v)->find_r_of_UV(u,v,M))
    # ϕ0 = Field(PS, (u,v)->0)
    # f0 = ((16*M^3)/r0)*exp(-r0/2M)
    # a0 = -sqrt(2*f0) 

    p  = 1e-1
    f0 = Field(PS, (u,v)->1)
    ϕ0 = Field(PS, (u,v)->p*exp(-(v-4)^2) + p*exp(-(u+7)^2))
    r0 = Field(PS, (u,v)->v-u)
    a0 = -sqrt(2*f0) 
    
    # Solve constraints at the incoming boundary and set boundary conditions
    bnda, bndr, bndϕ = initialdatasolver(a0, r0, ϕ0)

    # Compute constraints before solve
    E1, E2 = C(a0, r0, ϕ0)
    @show L2(E1), L2(E2)

    # Start solve
    fsolved, rsolved, ϕsolved = reshapeToTuple(PS, nlsolve(f!, reshapeFromTuple((a0, r0, ϕ0)); 
                                                           autodiff=:forward, show_trace=false, ftol=1e-9).zero)
    # Compute constraints after solve
    E1, E2 = C(fsolved, rsolved, ϕsolved)
    @show L2(E1), L2(E2)
   
    return (fsolved, rsolved, ϕsolved, E1, E2)
end


#--------------------------------------------------------------------
# Solve for an arbitrary spacetime and check for convergence 
# in L2 norm of constraints
#--------------------------------------------------------------------

NT = 10
L2CU = zeros(NT)
L2CV = zeros(NT)
NVEC = zeros(NT)

struct U end
struct V end

import DoubleFloats.Double64

for N in 4:NT
    @show N
    PS = ProductSpace(ChebyshevGL{U, N, Float64}(-8, -6), 
                      ChebyshevGL{V, N, Float64}( 3,  5))
    f, r, ϕ, cu, cv = nonlinearsolver(PS)
    L2CU[N] = L2(cu)
    L2CV[N] = L2(cv)
    NVEC[N] = N
end

using PyPlot
plot(NVEC, log10.(L2CU), "-o")
plot(NVEC, log10.(L2CV), "-o")
show()
