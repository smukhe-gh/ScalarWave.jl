#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Compactified coordinates Metric functions
#--------------------------------------------------------------------

function computericciscalar(grid::Grid, index::CartesianIndex)
    r       = grid.r_of_UV[index]
    t       = grid.t_of_UV[index]
    drdU    = grid.drdU[index]
    drdV    = grid.drdV[index]
    ddrdUdV = grid.ddrdUdV[index]
    M       = grid.params.mass
    riccis = (64*power(M,4) + 6*drdU*drdV*exp(r/(2.*M))*M*r + 
            6*ddrdUdV*exp(r/(2.*M))*M*power(r,2) - 
            ddrdUdV*exp(r/(2.*M))*power(r,3))/(32.*power(M,4)*power(r,2))
end

function computericci(grid::Grid, index::CartesianIndex)
    r       = grid.r_of_UV[index]
    t       = grid.t_of_UV[index]
    drdU    = grid.drdU[index]
    drdV    = grid.drdV[index]
    ddrdUdU = grid.ddrdUdU[index]
    ddrdUdV = grid.ddrdUdV[index]
    ddrdVdV = grid.ddrdVdV[index]
    M       = grid.params.mass
    theta   = grid.params.theta

    R11 = -((2*ddrdVdV*r + power(drdV,2)*(2 + r/M))/power(r,2))
    R22 = (16*power(M,3) + drdU*drdV*exp(r/(2.*M))*r + 
           ddrdUdV*exp(r/(2.*M))*power(r,2))/(16.*power(M,3))
    R33 = ((16*power(M,3) + drdU*drdV*exp(r/(2.*M))*r + 
            ddrdUdV*exp(r/(2.*M))*power(r,2))*power(sin(theta),2))/(16.*power(M,3))
    R41 = (-2*drdU*drdV*M + ddrdUdV*r*(-2*M + r))/(2.*M*power(r,2))
    R44 = -((2*ddrdUdU*r + power(drdU,2)*(2 + r/M))/power(r,2))

    return (R11, R22, R33, R41, R44)
end

function computedeterminant(grid::Grid, index::CartesianIndex)
    r = grid.r_of_UV[index]
    t = grid.t_of_UV[index]
    M = grid.params.mass 
    theta = grid.params.theta 
    determinant = (-1024*power(M,6)*power(r,2)*power(sin(theta),2))/exp(r/M)
end

function computemetric(grid::Grid, index::CartesianIndex)
    r = grid.r_of_UV[index]
    t = grid.t_of_UV[index]
    M = grid.params.mass 
    theta = grid.params.theta

    g01 = g10 = (-32*power(M,3))/(exp(r/(2.*M))*r)
    g22 = power(r,2) 
    g33 = power(r,2)*power(sin(theta),2)
    return (g01, g10, g22, g33)
end

function creategrid(params::Params)::Grid
    t_of_UV  = zeros(params.size)
    r_of_UV  = zeros(params.size)
    drdU     = zeros(params.size)
    drdV     = zeros(params.size)
    ddrdUdU  = zeros(params.size)
    ddrdUdV  = zeros(params.size)
    ddrdVdV  = zeros(params.size)
    UV       = createmesh(params) 

    Nx, Ny   = params.size
    @assert Nx == Ny
    N = Nx

    # compute (r,t) from (U,V)
    for index in CartesianRange(params.size)
        U,V = UV[index]
        t_of_UV[index], 
        r_of_UV[index] = find_TR_from_UV(U, V) 
    end

    # compute derivatives of (r,t) wrt (U,V)
    for index in CartesianRange(params.size)
        i = index[1]
        j = index[2]
        drdU[index]    = sum(chebd(i,m,N-1)*r_of_UV[m,j] for m in 1:N) 
        drdV[index]    = sum(chebd(j,n,N-1)*r_of_UV[i,n] for n in 1:N)
        ddrdUdU[index] = sum(chebd(i,l,N-1)*chebd(l,m,N-1)*r_of_UV[m,j] for m in 1:N, l in 1:N)
        ddrdUdV[index] = sum(chebd(j,n,N-1)*chebd(i,k,N-1)*r_of_UV[k,n] for k in 1:N, n in 1:N)
        ddrdVdV[index] = sum(chebd(j,k,N-1)*chebd(k,n,N-1)*r_of_UV[i,n] for n in 1:N, k in 1:N)
    end

    return Grid(params, t_of_UV, r_of_UV, 
                drdU, drdV, ddrdUdV, ddrdUdU, ddrdVdV)
end

function setmetric(grid::Grid)::Metric
    detg   = zeros(grid.params.size)
    riccis = zeros(grid.params.size)
    g01    = zeros(grid.params.size)
    g10    = zeros(grid.params.size)
    g22    = zeros(grid.params.size)
    g33    = zeros(grid.params.size)
    R11    = zeros(grid.params.size)
    R22    = zeros(grid.params.size)
    R33    = zeros(grid.params.size)
    R41    = zeros(grid.params.size)
    R44    = zeros(grid.params.size)

    # set all metric components and related quantities
    for index in CartesianRange(grid.params.size)
        detg[index] = computedeterminant(grid, index)
        riccis[index] = computericciscalar(grid, index)
        g01[index], g10[index],
        g22[index], g33[index] = computemetric(grid, index) 
        R11[index], R22[index], R33[index]
        R41[index], R44[index] = computericci(grid, index) 
    end

    return Metric(g01, g10, g22, g33, detg, riccis, R11, R22, R33, R41, R44)
end

function computeaction(grid::Grid, metric::Metric)::Float64
   riccis = metric.riccis
   detg   = metric.detg
   L      = riccis.*sqrt.(-detg)
   action = chebweights(grid.params.size[1]-1)'*L*chebweights(grid.params.size[2]-1)
end

