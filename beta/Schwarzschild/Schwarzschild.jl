#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Compactified coordinates Metric functions
#--------------------------------------------------------------------

function computericciscalar(grid::Grid, index::CartesianIndex)::Float64
    r       = grid.r[index]
    t       = grid.t[index]
    drdU    = grid.drdU[index]
    drdV    = grid.drdV[index]
    ddrdUdV = grid.ddrdUdV[index]
    M       = grid.params.mass
    scalar  = 64*power(M,4) + 6*drdU*drdV*exp(r/(2.*M))*M*r + 
              6*ddrdUdV*exp(r/(2.*M))*M*power(r,2) - 
              ddrdUdV*exp(r/(2.*M))*power(r,3)/(32.*power(M,4)*power(r,2))
    return scalar
end

function computemetricdet(grid::Grid, index::CartesianIndex)::Float64
    r = grid.r[index]
    t = grid.t[index]
    M = grid.params.mass 
    detg = (-1024*power(M,6)*power(r,2))/exp(r/M)
    return detg
end

function computemetric(grid::Grid, index::CartesianIndex)::Tuple
    r = grid.r[index]
    t = grid.t[index]
    M = grid.params.mass 
    g01 = g10 = (-32*power(M,3))/(exp(r/(2.*M))*r)
    g22 = g33 = power(r,2) 
    return (g01, g10, g22, g33)
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
        varlist.scalar[index] = computericciscalar(grid, index)
        varlist.g01[index], 
        varlist.g10[index],
        varlist.g22[index],
        varlist.g33[index]    = computemetric(grid, index) 
    end

    return varlist
end
