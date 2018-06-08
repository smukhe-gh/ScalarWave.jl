#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Compactified coordinates Metric functions
#--------------------------------------------------------------------

function setgrid(params::Params)::Grid
    grid = Grid(params,
                Array{Float64,2}(params.size),
                Array{Float64,2}(params.size),
                Array{Float64,2}(params.size),
                Array{Float64,2}(params.size),
                Array{Float64,2}(params.size),
                Array{Float64,2}(params.size),
                Array{Float64,2}(params.size),
                Array{Float64,2}(params.size),
                Array{Float64,2}(params.size),
                Array{Float64,2}(params.size),
                Array{Float64,2}(params.size),
                Array{Float64,2}(params.size),
                Array{Float64,2}(params.size),
                Array{Float64,2}(params.size)) 

    umin = params.xmin[1]
    umax = params.xmax[1]
    vmin = params.xmin[2]
    vmax = params.xmax[2]
    px   = params.p[1]
    py   = params.p[2]

    # create UV grid
    for index in CartesianRange(params.size)
        grid.U[index] = ((umax - umin)/2*chebx(index[1], px)) + (umax + umin)/2
        grid.V[index] = ((vmax - vmin)/2*chebx(index[2], py)) + (vmax + vmin)/2
    end

    # create TR grid
    for index in CartesianRange(params.size)
        grid.t[index], grid.r[index] = find_TR_of_UV(grid.U[index], grid.V[index]) 
    end

    drdU = dX_of_var(grid.r, grid, du) 
    drdV = dX_of_var(grid.r, grid, dv) 
    dtdU = dX_of_var(grid.t, grid, du) 
    dtdV = dX_of_var(grid.t, grid, dv) 
    ddrdUdU = ddX_of_var(grid.r, grid, du, du)
    ddrdUdV = ddX_of_var(grid.r, grid, du, dv)
    ddrdVdV = ddX_of_var(grid.r, grid, dv, dv)

    # it's time to compute the derivatives [FIX]
    for index in CartesianRange(params.size)
        grid.drdU[index] = drdU[index] 
        grid.drdV[index] = drdV[index]
        grid.dtdU[index] = dtdU[index]
        grid.dtdV[index] = dtdV[index]
        grid.ddrdUdU[index] = ddrdUdU[index] 
        grid.ddrdUdV[index] = ddrdUdV[index]
        grid.ddrdVdV[index] = ddrdVdV[index]
    end

    return grid
end

