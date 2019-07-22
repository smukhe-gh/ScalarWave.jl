#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 07-2019
# Test convergence with the operators
#--------------------------------------------------------------------

struct U end
struct V end

using IterativeSolvers
using Preconditioners

function computeawayfromaxis(N)
    SU  = ChebyshevGL{U, N, Float64}( 2, 4)
    SV  = ChebyshevGL{V, N, Float64}( 7, 9)
    SUV = ProductSpace(SU, SV)
    
    DU, DV = derivative(SUV)
    B = incomingboundary(SUV)
    I = identity(SUV)
    
    r = Field(SUV, (u,v)-> u-v)
    x = Field(SUV, (u,v)-> sinpi(u+v))
    y = Field(SUV, (u,v)->-(pi^2)*sinpi(u+v))
    
    P = DU*DV + ((DU*r)/r)*DV + ((DV*r)/r)*DU 
    return L2(P*x - y)
end

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
    P = enforceregularityonaxis(P, 3*DU*DV - DU*DU - DV*DV)
    return L2(P*x - y)
end

N  = []
EO = []
EA = []

using ProgressMeter
@showprogress for n in 1:30
    append!(N, n)
    append!(EA, computeawayfromaxis(n))
    append!(EO, computeonaxis(n))
end

using PyPlot
plot(N, log10.(EA), label="Away from axis")
plot(N, log10.(EO), label="On axis")
legend(frameon=false)
grid(true)
show()
