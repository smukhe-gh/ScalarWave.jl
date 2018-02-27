#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function prolongate(patch::Patch, M::Int)::Array{Float64,2}
    Nx = size(patch.value)[1] - 1
    Ny = size(patch.value)[2] - 1
    invndmx = inv(vandermonde(Nx))
    invndmy = inv(vandermonde(Ny))
    Px = zeros((Nx+1)*M, Nx+1)
    Py = zeros(Ny+1, (Ny+1)*M)
    for k in 1:M
        xg = chebgrid(Nx, M, k)
        yg = chebgrid(Ny, M, k)
        pvndmx = pseudovandermonde(Nx, xg) 
        pvndmy = pseudovandermonde(Ny, yg)
        Px[1+(k-1)*Nx:1+k*Nx, :] = pvndmx*invndmx 
        Py[:, 1+(k-1)*Ny:1+k*Ny] = pvndmy*invndmy 
    end
    return Px*patch.value*Py 
end

function restrict(dbase::Dict{Array{Int, 1}, Patch}, M::Int)::Array{Float64,2}
    Nx = size(dbase[[1,1]].value)[1] - 1
    Ny = size(dbase[[1,1]].value)[2] - 1
    invndmx = inv(vandermonde(Nx))
    invndmy = inv(vandermonde(Ny))
    Px = zeros((Nx+1)*M, Nx+1)
    Py = zeros(Ny+1, (Ny+1)*M)
    for k in 1:M
        xg = chebgrid(Nx, M, k)
        yg = chebgrid(Ny, M, k)
        pvndmx = pseudovandermonde(Nx, xg) 
        pvndmy = pseudovandermonde(Ny, yg) 
        Px[1+(k-1)*Nx:1+k*Nx, :] = pvndmx*invndmx 
        Py[:, 1+(k-1)*Ny:1+k*Ny] = pvndmy*invndmy 
    end
    Rx = pinv(Px)
    Ry = pinv(Py)
    gridval = zeros((Nx+1)*M, (Ny+1)*M) 
    for m in 1:M, n in 1:M
        gridval[1+(m-1)*Nx:1+m*Nx, 1+(n-1)*Ny:1+n*Ny] = dbase[[m,n]].value
    end
    return Rx*gridval*Ry
end


