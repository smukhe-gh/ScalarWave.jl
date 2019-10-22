#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2019
# Simulate Minkowski spacetime on axis
#--------------------------------------------------------------------

PS = ProductSpace(ChebyshevGL{U, 4, Float64}(0, 1), 
                  ChebyshevGL{V, 4, Float64}(0, 1))

A = axisboundary(PS)
B = incomingboundary(PS)
I = identity(PS)

r = Field(PS, (u,v)->2/(v-u))
w = Field(PS, (u,v)->4)
u = Field(PS, (u,v)->1)

display(mix(r, A, w))
@test mix(r,A,w) ≈ mixonaxis(r, w)

display(mix(w, B, u))
@test mix(w, B, u) ≈ mixonincomingboundary(w, u)
