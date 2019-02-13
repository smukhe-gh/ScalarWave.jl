#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# Plot fields using PyPlot
#--------------------------------------------------------------------

uu = Field(GaussLobatto{U, 10, 3.0, -3.0}, x->sin(x))
vv = Field(ProductSpace{GaussLobatto{V, 8, 3, -3},
                        GaussLobatto{U, 6, 5, -5}}, (x,y)->exp(-x^2) + exp(-y^2)) 
nest = refine(uu)
nest[[2]] = refine(nest[[2]])
nest[[2]][[2]] = refine(nest[[2]][[2]])

nest2D = refine(vv)
nest2D[[2,2]] = refine(nest2D[[2,2]])

nest2Dspace = refine(ProductSpace{GaussLobatto{V, 8, 3, -3},
                                  GaussLobatto{U, 6, 5, -5}})
nest2Dspace[[2,2]] = refine(nest2Dspace[[2,2]])

@testset "PyPlot" begin
    @test plot(uu) == 0
    @test plot(uu, 101) == 0
    @test contour(vv, 101) == 0
    @test contourf(vv, 101) == 0
    @test pcolormesh(vv) == 0
    @test pcolormesh(vv, 101) == 0
    @test plot(nest, 100) == 0
    @test contour(nest2D, 100, globalmax=maximum(nest2D), globalmin=minimum(nest2D), globallevels=levels(nest2D, globallength=20)) == 0
    @test contourf(nest2D, 100, globalmax=maximum(nest2D), globalmin=minimum(nest2D), globallevels=levels(nest2D, globallength=20)) == 0
    close()
end

contourf(nest2D, 100, globalmax=maximum(nest2D), globalmin=minimum(nest2D), globallevels=levels(nest2D, globallength=20))
plot(nest2Dspace)
show()
