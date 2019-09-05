#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 07-2019
# Test functions inside functions
#--------------------------------------------------------------------

using NLsolve

function nonlinearsolver(PS::ProductSpace{S1, S2})::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}

    # Compute constraints
    function C(f::Field{S}, r::Field{S}, ϕ::Field{S})::NTuple{2, Field{S}} where {S}
        E1 = 2*(DU*DU*r - (1/f)*(DU*f)*(DU*r)) + r*(DU*ϕ)^2
        E2 = 2*(DV*DV*r - (1/f)*(DV*f)*(DV*r)) + r*(DV*ϕ)^2
        return (E1, E2)
    end

    # Compute residuals
    function F(f::Field{S}, r::Field{S}, ϕ::Field{S})::NTuple{3, Field{S}} where {S}
        resf = DU*DV*log(abs(f)) + (2/r)*(DU*DV*r) + 2*(DU*ϕ)*(DV*ϕ)
        resr = 2*DU*(DV*r) +  (2/r)*(DU*r)*(DV*r) + (f/r)
        resϕ = DU*DV*ϕ + (1/r)*(DU*r)*(DV*ϕ) + (1/r)*(DV*r)*(DU*ϕ)
    
        resbndf = B*(f - bndf)  
        resbndr = B*(r - bndr)  
        resbndϕ = B*(ϕ - bndϕ)  
    
        finalresf = (I-B)*resf + resbndf
        finalresr = (I-B)*resr + resbndr
        finalresϕ = (I-B)*resϕ + resbndϕ
    
        return (finalresf, finalresr, finalresϕ)
    end
    
    # Wrapper for nlsolve
    function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
        fvec[:] = reshapeFromTuple(F(reshapeToTuple(PS, x)...))
    end

    # Compute common operators
    DU, DV = derivative(PS)
    B = incomingboundary(PS)
    I = identity(PS)

    # Set up free variables
    ϕ0 = Field(PS, (u,v)->exp(-(v-4)^2))  # Incoming wave travelling in the u-direction
    r0 = Field(PS, (u,v)->v-u)
    f0 = Field(PS, (u,v)->1)

    # Schwarzschild spacetime 
    M = 1.0
    r0 = Field(PS, (u,v)->find_r_of_UV(u,v,M))
    f0 = ((16*M^3)/r0)*exp(-r0/2M)
    ϕ0 = Field(PS, (u,v)->0)
    
    # Solve constraints at the incoming boundary and set boundary conditions
    f0, r0, ϕ0 = initialdatasolver(f0, r0, ϕ0)
    bndf = B*f0
    bndr = B*r0
    bndϕ = B*ϕ0

    # Compute constraints before solve
    E1, E2 = C(f0, r0, ϕ0)
    @show L2(E1), L2(E2)

    # Start solve
    fsolved, rsolved, ϕsolved = reshapeToTuple(PS, nlsolve(f!, reshapeFromTuple((f0, r0, ϕ0)); 
                                                           autodiff=:forward, show_trace=true, ftol=1e-9).zero)
    # Compute constraints after solve
    E1, E2 = C(fsolved, rsolved, ϕsolved)
    @show L2(E1), L2(E2)
   
    return (fsolved, rsolved, ϕsolved)
end


#--------------------------------------------------------------------
# Solve for an arbitrary spacetime
#--------------------------------------------------------------------

struct U end
struct V end
S1 = ChebyshevGL{U, 23, Float64}(-8, -6)
S2 = ChebyshevGL{V, 23, Float64}( 3,  5)
PS = ProductSpace(S1, S2)

f, r, ϕ = nonlinearsolver(PS)
