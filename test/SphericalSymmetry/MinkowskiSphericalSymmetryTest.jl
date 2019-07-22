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

using IterativeSolvers
using Preconditioners

SU  = ChebyshevGL{U, 20, Float64}(-1, 1)
SV  = ChebyshevGL{V, 20, Float64}(-1, 1)
SUV = ProductSpace(SU, SV)

u = Field(SUV, (u,v)->u)
v = Field(SUV, (u,v)->v)
r = Field(SUV, (u,v)->u-v)

DU, DV = derivative(SUV)
B = incomingboundary(SUV)
I = identity(SUV)

u = Field(SUV, (u,v)-> sin(u+v))
b = Field(SUV, (u,v)->-sin(u+v))
P = DU*DV + ((DU*r)/r)*DV + ((DV*r)/r)*DU 

if false
    L = enforceregularityonaxis(P, DU*DV)
    A = reshape(L ⊕ B)
    c = reshape(B*u + (I-B)*b)
    p = DiagonalPreconditioner(A) 

    # ϕ = reshape(SUV, gmres(A,c; verbose=true))
    ϕ = reshape(SUV, A \ c)
    # ϕ = reshape(SUV, cg(A,c, Pl=p))

    @show L1(L*ϕ - b)
    @show L1(ϕ - u)
    pcolormesh(ϕ - u)
    show()
end

# Let's assume the limit exists and replace the singular 
# part of the operator with a regular version. 

if true
    Q = enforceregularityonaxis(P, 3*DU*DV - DU*DU - DV*DV)
    # Q = enforceregularityonaxis(P, DU*DV)
    @show cond(Q ⊕ B)
    ϕ = solve(Q ⊕ B, B*u + (I-B)*b)

    @show L1(Q*u - b)
    @show L1(Q*ϕ - b)
    @show L1(ϕ - u)

    pcolormesh(ϕ)
    show()
end
