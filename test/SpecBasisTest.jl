#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
#--------------------------------------------------------------------

S = GaussLobatto(U, 4, 12, -2)
x = Field(S, u->u^3)
y = Field(S, u->3*u^2)
D = derivative(S)

@test D*x ≈ y 

