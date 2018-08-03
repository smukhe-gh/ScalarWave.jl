#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 05-2018
#--------------------------------------------------------------------

parameters = Dict("px"   => 20,
                  "py"   => 20,
                  "umin" => -20.0,
                  "umax" => -0.1,
                  "vmin" => 0.1,
                  "vmax" => 20.0,
                  "mass" => 1.0)

params  = dict2struct(parameters)
grid    = setgrid(params)
varlist = setvarlist(grid) 
dXdU    = 2/(grid.params.dmax[1] - grid.params.dmin[1])
dXdV    = 2/(grid.params.dmax[2] - grid.params.dmin[2])

# test coordinates
for index in CartesianRange(params.size)
    U = grid.U[index]
    V = grid.V[index]
    t = grid.t[index]
    r = grid.r[index]
    @test U ≈ find_UV_of_TR(t,r, grid.params.mass)[1] 
    @test V ≈ find_UV_of_TR(t,r, grid.params.mass)[2] 
end

# test incoming scalarwave
dict = SchwarzschildDistribute(u->0, v->exp(-v^2/0.1), (u,v)->0, parameters)

# plot the solution
#drawmultipatch(dict, "schwarzschild-test-wave-operator")
