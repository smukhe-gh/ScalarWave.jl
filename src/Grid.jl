#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function restrictmodes!(coeffs::Array{Float64,1}, M::Int)::Array{Float64,1}
    if M+1 > size(coeffs)[1]
        error("Cannot restrict to larger number of modes")
    else
        for index in CartesianRange(size(coeffs))
            index[1] > M+1 ? coeffs[index] = 0 : coeffs[index] = coeffs[index] 
        end
    end
    return coeffs
end

function restrictmodes!(coeffs::Array{Float64,2}, M::Int, N::Int)::Array{Float64,2}
    if M+1 > size(coeffs)[1] || N+1 > size(coeffs)[2]
        error("Cannot restrict to larger number of modes")
    else    
        for index in CartesianRange(size(coeffs))
            index[1] > M + 1 || index[2] > N + 1 ? coeffs[index] = 0 : coeffs[index] = coeffs[index] 
        end
    end
    return coeffs
end

function prolongatemodes(coeffs::Array{Float64,1}, M::Int)::Array{Float64,1}
    if M+1 < size(coeffs)[1]
        error("Cannot prolongate to smaller  number of modes")
    else
        ncoeffs = zeros(M+1)
        oldM = size(coeffs)[1] - 1
        for index in CartesianRange(size(ncoeffs))
            index[1] < oldM+1 ? ncoeffs[index] = coeffs[index] : ncoeffs[index] = ncoeffs[index] 
        end
    end
    return ncoeffs
end

function prolongatemodes(coeffs::Array{Float64,2}, M::Int, N::Int)::Array{Float64,2}
    if M+1 < size(coeffs)[1] || N+1 < size(coeffs)[2] 
        error("Cannot prolongate to smaller number of modes")
    else
        ncoeffs = zeros(M+1,N+1)
        oldM = size(coeffs)[1] - 1
        oldN = size(coeffs)[2] - 1
        for index in CartesianRange(size(ncoeffs))
            index[1] < oldM+1 && index[2] < oldN+1 ? ncoeffs[index] = coeffs[index] : ncoeffs[index] = ncoeffs[index] 
        end
    end
    return ncoeffs
end

function prolongateOP(Nx::Int, M::Int)::Array{Float64,2}
    invndmx = inv(vandermonde(Nx))
    Px = zeros((Nx+1)*M, Nx+1)
    for k in 1:M
        pvndmx = pseudovandermonde(Nx, chebgrid(Nx, M, k))
        Px[1+(k-1)*(Nx+1):k*(Nx+1), :] = pvndmx*invndmx
    end
    return Px
end 

function restrictOP(Nx::Int, M::Int)::Array{Float64,2}
    return pinv(prolongateOP(Nx, M))
end 

function prolongatePatch(patch::Patch, M::Int)::Dict{Array{Int, 1}, Patch}
    Nx = size(patch.value)[1] - 1
    Ny = size(patch.value)[2] - 1
    Px = prolongateOP(Nx, M) 
    Py = prolongateOP(Ny, M)
    mPatch = Px*patch.value*Py'
    return array2dict(mPatch, Nx, Ny, M)
end

function restrictPatch(dbase::Dict{Array{Int, 1}, Patch})::Patch
    Nx = size(dbase[[1,1]].value)[1] - 1
    Ny = size(dbase[[1,1]].value)[2] - 1
    M  = convert(Int, sqrt(length(dbase))) 
    Rx = restrictOP(Nx, M)
    Ry = restrictOP(Ny, M)
    sPatch = dict2array(dbase) 
    return Patch([1,1], Rx*sPatch*Ry')
end

#--------------------------------------------------------------------
# Additional utility functions for computing coordinates
# and their derivatives on the grid
#--------------------------------------------------------------------

function find_TR_of_UV(U::Float64, V::Float64)::Tuple
    @vars x
    u = tan(find_zero(atan((tan(x)/(sqrt(1+tan(x)^2)))*log(1+tan(x)^2)) - U, (-pi/2, pi/2)))
    v = tan(find_zero(atan((tan(x)/(sqrt(1+tan(x)^2)))*log(1+tan(x)^2)) - V, (-pi/2, pi/2)))
    t = u+v
    r = v-u
    return (t, r)
end

function find_UV_of_TR(t::Float64, r::Float64)::Tuple
    u = (t-r)/2
    v = (t+r)/2
    U = atan((u/(sqrt(1+u^2)))*log(1+u^2))
    V = atan((v/(sqrt(1+v^2)))*log(1+v^2))
    return (U,V)
end

function dX_of_var(var::Array{Float64,2}, grid::Grid, X::D)
    dvar = zeros(grid.params.size)
    dXdU = 2/(grid.params.dmax[Int(X)] - grid.params.dmin[Int(X)])
    for index in CartesianRange(grid.params.size)
        (i, j) = index.I
        if Int(X) == 1
            dvar[index] = dXdU*sum(chebd(i, m, grid.params.mode[Int(X)])*var[m, j] for m in 1:grid.params.size[Int(X)])
        else
            dvar[index] = dXdU*sum(chebd(j, n, grid.params.mode[Int(X)])*var[i, n] for n in 1:grid.params.size[Int(X)])
        end
    end
    return dvar
end

function ddX_of_var(var::Array{Float64,2}, grid::Grid, X1::D, X2::D)
    ddvar = zeros(grid.params.size)
    dXdU  = 2/(grid.params.dmax[Int(X1)] - grid.params.dmin[Int(X1)])
    dXdV  = 2/(grid.params.dmax[Int(X2)] - grid.params.dmin[Int(X2)])
    for index in CartesianRange(grid.params.size)
        i, j = index[1], index[2]
        if (Int(X1) == 1 && Int(X2) == 1)
            ddvar[index] = dXdU*dXdU*sum(chebd(i, l, grid.params.mode[Int(X1)])*chebd(l, m, grid.params.mode[Int(X1)])*var[m,j] 
                                         for m in 1:grid.params.size[Int(X1)], l in 1:grid.params.size[Int(X1)])
        elseif (Int(X1) == 2 && Int(X2) == 2)
            ddvar[index] = dXdV*dXdV*sum(chebd(j, k, grid.params.mode[Int(X2)])*chebd(k, n, grid.params.mode[Int(X2)])*var[i,n] 
                                         for n in 1:grid.params.size[Int(X2)], k in 1:grid.params.size[Int(X2)])
        else
            ddvar[index] = dXdU*dXdV*sum(chebd(j, n, grid.params.mode[Int(X1)])*chebd(i, k, grid.params.mode[Int(X2)])*var[k,n] 
                                         for k in 1:grid.params.size[Int(X1)], n in 1:grid.params.size[Int(X2)])
        end
    end
    return ddvar
end
