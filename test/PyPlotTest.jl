#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# Plot fields using PyPlot
#--------------------------------------------------------------------

# test 1D fields
uu = Field(GaussLobatto{U, 10, 3.0, -3.0}, x->sin(x))
plot(uu)
plot(uu, collocation=true)
close()

# test 2D fields
vv = Field(ProductSpace{GaussLobatto{V, 8, 3, -1},
                        GaussLobatto{U, 6, 5, -1}}, (x,y)->(x+y)) 
contourf(vv)
show()
