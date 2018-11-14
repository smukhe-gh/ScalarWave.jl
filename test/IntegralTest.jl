#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Test integration on a patch
#--------------------------------------------------------------------

# Test 1D operation
P = 80
S = GaussLobatto(U, P, 3*pi, 4)
F = Field(S, U->sin(U)^2)
W = integral(S)  

@show W*F

# test 2D operations

P1, P2 = 28, 26
SUV = ProductSpace{GaussLobatto(V, P2, 8*pi, 4), GaussLobatto(U, P1, -3, -5*pi)}
FUV = Field(SUV, (U,V)->U^2 + V^4)
WUV = integral(SUV)
@show WUV*FUV

