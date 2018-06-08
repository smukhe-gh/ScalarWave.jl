#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Compactified coordinates Metric functions
#--------------------------------------------------------------------

function omega(grid::Grid, index::CartesianIndex)::Float64
    r = grid.r[index]
    t = grid.t[index]
    u = (t - r)/2
    v = (t + r)/2
    return exp((t-r)/2M)*((2*power(u,2) + Log(1 + power(u,2)))*(2*power(v,2) + Log(1 + 
            power(v,2))))/(sqrt(1 + power(u,2))*sqrt(1 + power(v,2))*(1 + power(u,2) + 
            power(u,2)*power(Log(1 + power(u,2)),2))*(1 + power(v,2) + 
            power(v,2)*power(Log(1 + power(v,2)),2)))
end

function SchwarzschildOP(grid::Grid, varlist::VarList)::Float64
   """
   XXX: Test for large values 
   """
end

function computericciscalar(grid::Grid, index::CartesianIndex)::Float64
    r       = grid.r[index]
    t       = grid.t[index]
    drdU    = grid.drdU[index]
    drdV    = grid.drdV[index]
    ddrdUdV = grid.ddrdUdV[index]
    M       = grid.params.mass
    # FIXME: The exp(r/M) is causing issues
    scalar  = (32 + (exp(r/M)*r*(3*drdU*drdV*M + ddrdUdV*(3*M 
                  - r)*r))/power(M,4))/(16.*power(r,2))
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
