#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 10-2019
# Test basis change functions
#--------------------------------------------------------------------

PS = ProductSpace(ChebyshevGL{U, 4, Float64}(0.0,  1.0), 
                  ChebyshevGL{V, 4, Float64}(0.0,  1.0))

I = identity(PS)
B = incomingboundary(PS)
u = Field(PS.S1, u->1+u^3)
x = Field(PS.S1, u->u)

ω = basistransform(u)
η = basistransform(basistransform(u))
@test η ≈ u

@show x
@show basistransform(u)

@show (0.8)^3 +  1
@show u(0.8)
