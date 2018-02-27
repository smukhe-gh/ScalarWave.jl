#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

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


