#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Compactified coordinates Metric functions
#--------------------------------------------------------------------

function computemetricdet(grid::Grid, index::CartesianIndex)::Float64
    r = grid.r[index]
    t = grid.t[index]
    M = grid.params.mass 
    detg = (-1024*power(M,6)*power(r,2))/exp((2*r)/M) 
    return detg
end

function computemetric(grid::Grid, index::CartesianIndex)::Tuple
    r = grid.r[index]
    t = grid.t[index]
    M = grid.params.mass 
    g01 = g10 = (-32*power(M,3))/(exp(r/M)*r) 
    g22 = g33 = power(r,2) 
    return (g01, g10, g22, g33)
end

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

    umin = params.dmin[1]
    umax = params.dmax[1]
    vmin = params.dmin[2]
    vmax = params.dmax[2]
    px   = params.mode[1]
    py   = params.mode[2]

    # create UV grid
    for index in CartesianRange(params.size)
        grid.U[index] = ((umax - umin)/2*chebx(index[1], px)) + (umax + umin)/2
        grid.V[index] = ((vmax - vmin)/2*chebx(index[2], py)) + (vmax + vmin)/2
    end

    # create TR grid
    for index in CartesianRange(params.size)
        grid.t[index], grid.r[index] = find_TR_of_UV(grid.U[index], grid.V[index], grid.params.mass) 
    end

    drdU = dX_of_var(grid.r, grid, D(u)) 
    drdV = dX_of_var(grid.r, grid, D(v)) 
    dtdU = dX_of_var(grid.t, grid, D(u)) 
    dtdV = dX_of_var(grid.t, grid, D(v)) 
    ddrdUdU = ddX_of_var(grid.r, grid, D(u), D(u))
    ddrdUdV = ddX_of_var(grid.r, grid, D(u), D(v))
    ddrdVdV = ddX_of_var(grid.r, grid, D(v), D(v))

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

function setvarlist(grid::Grid)::VarList
    varlist = VarList(zeros(grid.params.size),
                      zeros(grid.params.size),
                      zeros(grid.params.size),
                      zeros(grid.params.size),
                      zeros(grid.params.size),
                      zeros(grid.params.size)) 

    for index in CartesianRange(grid.params.size)
        varlist.detg[index]   = computemetricdet(grid, index)
        varlist.g01[index], 
        varlist.g10[index],
        varlist.g22[index],
        varlist.g33[index]    = computemetric(grid, index) 
    end

    return varlist
end

function SchwarzschildDerivOP(grid::Grid, varlist::VarList)::Array{Float64, 4}
    Nx = grid.params.mode[1]
    Ny = grid.params.mode[2]
    @assert Nx == Ny
    dXdU = 2/(grid.params.dmax[1] - grid.params.dmin[1])
    dXdV = 2/(grid.params.dmax[2] - grid.params.dmin[2])
	operator = zeros(Nx+1, Ny+1, Nx+1, Ny+1)
    detg     = varlist.detg
    g01      = varlist.g01
    for index in CartesianRange(size(operator))    
        k = index.I[1]
        i = index.I[2]
        l = index.I[3]
        j = index.I[4]   
        operator[index] = (dXdU*dXdV)^2*(2*g01[i,k]*chebw(i,Ny)*chebw(k,Nx)*chebd(k,l,Nx)*chebd(i,j,Ny))*sqrt(-detg[i,k])	
	end
	return operator
end
