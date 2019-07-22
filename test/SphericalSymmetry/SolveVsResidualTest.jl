#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 07-2019
# Solve the wave equation with a fixed background 
# in spherically-symmetric spacetimes 
#--------------------------------------------------------------------

struct U end
struct V end

function computeonaxis(N)
    SU  = ChebyshevGL{U, N, Float64}(-1, 1)
    SV  = ChebyshevGL{V, N, Float64}(-1, 1)
    SUV = ProductSpace(SU, SV)
    
    DU, DV = derivative(SUV)
    B = incomingboundary(SUV)
    I = identity(SUV)
    
    r = Field(SUV, (u,v)-> u-v)
    x = Field(SUV, (u,v)-> sinpi(u+v))
    y = Field(SUV, (u,v)->-(pi^2)*sinpi(u+v))
    
    P = DU*DV + ((DU*r)/r)*DV + ((DV*r)/r)*DU 
    Q = enforceregularityonaxis(P, 3*DU*DV - DU*DU - DV*DV)
    ϕ = solve(Q ⊕ B, B*x + (I-B)*y)
    return (L2(Q*x-y), L2(ϕ - x))
end

N  = []
AE = []
NE = []

using ProgressMeter
@showprogress for n in 1:30
    append!(N, n)
    ae, ne = computeonaxis(n)
    append!(AE, ae)
    append!(NE, ne)
end

using PyPlot
plot(N, log10.(NE), label="solve")
plot(N, log10.(AE), label="exact residual")
legend(frameon=false)
grid(true)
show()
