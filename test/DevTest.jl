#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 07-2019
# Solve the wave equation with a fixed background 
# in spherically-symmetric spacetimes 
#--------------------------------------------------------------------

struct U end
struct V end

#--------------------------------------------------------------------
# Test away from the origin
#--------------------------------------------------------------------

if false
    SU  = ChebyshevGL{U, 20, Float64}(-1, 1)
    SV  = ChebyshevGL{V, 20, Float64}( 2, 4)
    SUV = ProductSpace(SU, SV)
    
    u = Field(SUV, (u,v)->u)
    v = Field(SUV, (u,v)->v)
    r = Field(SUV, (u,v)->u-v)
    
    DU, DV = derivative(SUV)
    B = incomingboundary(SUV)
    I = identity(SUV)
    
    L = DU*DV + ((DU*r)/r)*DV + ((DV*r)/r)*DU 
    
    u = Field(SUV, (u,v)->sin(u+v))
    b = Field(SUV, (u,v)->sin(u+v))
    @show L1(L*u + b)
    
    ϕ  = solve(L ⊕ B, B*u - (I-B)*b)
    @show L1(L*ϕ + b)
end


#--------------------------------------------------------------------
# Test on axis
#--------------------------------------------------------------------

if false
    SU  = ChebyshevGL{U, 30, Float64}(-1, 1)
    SV  = ChebyshevGL{V, 30, Float64}(-1, 1)
    SUV = ProductSpace(SU, SV)
    
    u = Field(SUV, (u,v)->u)
    v = Field(SUV, (u,v)->v)
    r = Field(SUV, (u,v)->u-v)
    
    DU, DV = derivative(SUV)
    B = incomingboundary(SUV)
    I = identity(SUV)
    
    L = DU*DV + ((DU*r)/r)*DV + ((DV*r)/r)*DU 
    L = enforceregularityonaxis(L, DU*DV)
    @show cond(L ⊕ B)
    
    u = Field(SUV, (u,v)->sin(u+v))
    b = Field(SUV, (u,v)->sin(u+v))
    @show L1(L*u + b)
    
    using IterativeSolvers
    using Preconditioners

    A = reshape(L ⊕ B)
    c = reshape(B*u - (I-B)*b)
    P = DiagonalPreconditioner(A) 

    # ϕ = reshape(SUV, gmres(A,c; verbose=true))
    # ϕ = reshape(SUV, A \ c)
    ϕ = reshape(SUV, cg(A,c, Pl=P))

    @show L1(L*ϕ + b)
    @show L1(ϕ - u)

    pcolormesh(ϕ - u)
    show()
end

