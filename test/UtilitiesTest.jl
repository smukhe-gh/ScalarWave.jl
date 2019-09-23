#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 10-2019
# Test utilities
#--------------------------------------------------------------------

PS = ProductSpace(ChebyshevGL{U, 20, Float64}(-1.0,  1.0), 
                  ChebyshevGL{V, 20, Float64}( 0.2,  0.9))

(a0, r0, ϕ0) = background(PS, (u,v)->exp(-v^2))
(abnd, rbnd, ϕbnd) = combineUVboundary.(computeUboundary((a0, r0, ϕ0)), computeVboundary((a0, r0, ϕ0)), :incoming)

@show lineconstraint(extractUboundary.((abnd, rbnd, ϕbnd), :incoming)...)
@show lineconstraint(extractVboundary.((abnd, rbnd, ϕbnd), :incoming)...)


# Schwarzschild spacetime 
M  = T(1.0)
r0 = Field(PS, (u,v)->find_r_of_UV(u,v,M))
ϕ0 = Field(PS, (u,v)->0)
f0 = ((16*M^3)/r0)*exp(-r0/2M)
a0 = -sqrt(2*f0) 

@show L2.(constraints(a0, r0, ϕ0))
