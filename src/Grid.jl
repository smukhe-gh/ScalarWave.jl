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
            index[1] > M+1 || index[2] > N+1 ? coeffs[index] = 0 : coeffs[index] = coeffs[index] 
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
