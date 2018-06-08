#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Compactified coordinates Metric functions
#--------------------------------------------------------------------

using PyPlot, ScalarWave

parameters = Dict("px"   => 80,   "py"   => 80,
                  "umin" => -1.5, "umax" => 1.5,
                  "vmin" => -1.5, "vmax" => 1.5,
                  "mass" => 1.0)

params  = dict2struct(parameters)
grid    = setgrid(params)

levels =[2, 5, 10, 20]
CS = contour(grid.U, grid.V, grid.r, colors="b")
clabel(CS, inline=1, fontsize=10)
CQ = contour(grid.V, grid.U, grid.t, colors="r")
clabel(CQ, inline=1, fontsize=10)
ylabel("U")
xlabel("V")
show()
