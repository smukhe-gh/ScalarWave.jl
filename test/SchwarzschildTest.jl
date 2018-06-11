#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 05-2018
#--------------------------------------------------------------------

parameters = Dict("px"   => 4,
                  "py"   => 4,
                  "umin" => 0.0,
                  "umax" => 0.5,
                  "vmin" => 0.2,
                  "vmax" => 1.7,
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
    @test U ≈ find_UV_of_TR(t,r)[1] 
    @test V ≈ find_UV_of_TR(t,r)[2] 
end

