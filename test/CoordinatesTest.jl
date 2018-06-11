#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 05-2018
#--------------------------------------------------------------------

U, V, M  = (-0.00001, 0.9, 2.4)
T, R     = find_TR_of_UV(U, V, M)
RU, RV   = find_UV_of_TR(T, R, M) 
@test RU ≈ U  
@test RV ≈ V  

parameters = Dict("px"   => 15,
                  "py"   => 15,
                  "umin" => -15.0,
                  "umax" => -0.1,
                  "vmin" => 0.1,
                  "vmax" => 15.0,
                  "mass" => 1.0)

params  = dict2struct(parameters)
grid    = setgrid(params)
varlist = setvarlist(grid) 

# test coordinates
for index in CartesianRange(params.size)
    U = grid.U[index]
    V = grid.V[index]
    t = grid.t[index]
    r = grid.r[index]
    @test U ≈ find_UV_of_TR(t,r, grid.params.mass)[1] 
    @test V ≈ find_UV_of_TR(t,r, grid.params.mass)[2] 
end

# show the grid you're working with.
# plotcoordgrid(params)
