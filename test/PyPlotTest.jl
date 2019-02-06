#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# Plot fields using PyPlot
#--------------------------------------------------------------------

uu = Field(GaussLobatto{U, 10, 3.0, -3.0}, x->sin(x))
vv = Field(ProductSpace{GaussLobatto{V, 8, 3, -1},
                        GaussLobatto{U, 6, 5, -1}}, (x,y)->(x+y)) 

plot(uu)
plot(uu, 100)

contour(vv, 100)
contourf(vv, 100)

pcolormesh(vv)
pcolormesh(vv, 100)

close()


