#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Test metric functions
#--------------------------------------------------------------------

function testexpansion(grid::Grid, metric::Metric)::Float64
    r = grid.r_of_UV[index]
    M = grid.params.mass 
    return -3/2\sqrt(2M/r^3)
end

params = Params((4,4),              # no. of points in U, V 
                (0.5, 0.8),       # span in U
                (1.3, 1.45),       # span in V
                1.0,                # BH mass
                pi/2)               # theta 

grid   = creategrid(params)
metric = setmetric(grid)

@show metric.riccis
@show metric.R11
@show metric.R22
@show metric.R33
@show metric.R41
@show metric.R44
@show computeaction(grid, metric)
