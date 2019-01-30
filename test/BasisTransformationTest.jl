#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# Basis transformation Test
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# test 1D spaces
#--------------------------------------------------------------------

f = Field(GaussLobatto{U, 10, 2, -3}, rand(11))
fbar = basistransform(basistransform(f))
@test f ≈ fbar
u = Field(GaussLobatto{U, 10, 2, -3}, x->x^3)

#--------------------------------------------------------------------
# test 2D spaces
#--------------------------------------------------------------------

S = ProductSpace{GaussLobatto{V, 9, 1, -1}, 
                 GaussLobatto{U, 4, 1, -1}}

uu = Field(S, (U,V)->U)
vv = Field(S, (U,V)->V)
ff = Field(S, (U,V)->U^2 + V^3)

cc  = basistransform(ff)
ffr = basistransform(cc) 
ccr = basistransform(ffr) 
@test ffr.value ≈ ff.value
@test ccr.value ≈ cc.value

W = Field(S, (U,V)->exp(U) + V)
@show W(0.3, 0.4)

S = ProductSpace{GaussLobatto{V, 29, 5, -2}, 
                 GaussLobatto{U, 34, 9, 6}}
W = Field(S, (U,V)->exp(U) + V)
@show W(6.5, -1)

S = ProductSpace{GaussLobatto{V, 29, 5, -2}, 
                 GaussLobatto{U, 34, 9, -6}}
W = Field(S, (U,V)->exp(U) + V)
@show W(0.3, 0.4)
