#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2019
# Test new functions for 1D spaces 
#--------------------------------------------------------------------

S = GaussLobatto{U, 3, 1, -1}
@test maximum(S) ==  1
@test minimum(S) == -1
@test order(S)   ==  3
@test length(S)  ==  4

x = Field(S, x->x)
y = Field(S, x->1)

@test derivative(S)*x ≈ y
@test tr(integral(S)*x).value ≈  0.0
@test identity(S)*x ≈ x 
@test (incomingboundary(S)*x).value[1] == x.value[1]
@test (outgoingboundary(S) x).value[end] == x.value[end]

