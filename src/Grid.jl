#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function prolongate(patch::Patch, M::Int)::Dict{Array{Float64,1}, Patch}
    Nx = size(patch.value)[1] - 1
    Ny = size(patch.value)[2] - 1
    xp = chebgrid(Nx) 
    yp = chebgrid(Ny)
    fmodal = extractPatchCoeffs(patch)  
    dbase  = Dict()
    for k in 1:M, l in 1:M
        loc = [k,l] 
        xg  = chebgrid(Nx, M, loc[1])
        yg  = chebgrid(Ny, M, loc[2]) 
        # for testing. we do not need to store these matrices.
        pvndmx = pseudovandermonde(Nx, xg) 
        pvndmy = pseudovandermonde(Ny, yg)
        for index in CartesianRange(size(fnodal))
            i = index.I[1]
            j = index.I[2]
            fnodal[i,j] = sum(pvndmx[i,m]*fmodal[m,n]*pvndmy[j,n] for m in 1:Nx+1, n in 1:Ny+1)
        end
        dbase[loc] = Patch(loc, fpatch) 
    end
    return dbase
end

function restrict(dbase::Dict{Array{Int, 1}, Patch}, M::Int)::Patch
    Nx = size(dbase[[1,1]])[1] - 1
    Ny = size(dbase[[1,1]])[2] - 1
    if Nx % M != 0 || Nx % M != 0
        error("Nx and Ny needs to be a multiple of M at the moment.")
    end
    dbase = Dict() 
    for i in 1:M, j in 1:M
        loc = [i,j]    
        fmodal = extractPatchCoeffs(dbase[[i,j]])
        xg = chebgrid(Nx/M, M, loc[1])
        yg = chebgrid(Ny/M, M, loc[2]) 
        pvndmx = pseudovandermonde(Nx, xg) 
        pvndmy = pseudovandermonde(Ny, yg)
        for index in CartesianRange(size(fnodal))
            i = index.I[1]
            j = index.I[2]
            fnodal[i,j] = sum(pvndmx[i,m]*fmodal[m,n]*pvndmy[j,n] for m in 1:Nx+1, n in 1:Ny+1)
        end
        dbase[loc] = Patch(loc, fpatch) 
    end
    
    fglobalgrid = pinv(cPx)*fgrid*pinv(cPy)
    return Patch([1,1], fglobalgrid) 
end


