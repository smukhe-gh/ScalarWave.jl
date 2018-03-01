#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function getPatchIC{T<:Integer}(fn::Function, s::T, Nx::T, Ny::T, M::T, loc::Array{T,1})::Boundary
    if s == 0
        return Boundary(0, projectonPatchBndbyRestriction(fn, Nx, M, loc[1]))
    elseif s == 1
        return Boundary(1, projectonPatchBndbyRestriction(fn, Ny, M, loc[2]))
    else
        error("Unknown direction passed.")
    end
end

function getPatchBnd(patch::Patch, s::Int)::Boundary
    if s == 0
        boundary = Boundary(0, patch.value[:, end])
    elseif s == 1
        boundary = Boundary(1, patch.value[end, :])
    else 
        error("Unknown direction passed")
    end
    return boundary
end

function calcPatch(bndx::Boundary, bndy::Boundary, RHS::Array{Float64,2}, 
                   derivOP::Array{Float64,4}, boundaryOP::Array{Float64,4},
                   loc::Array{Int,1})::Patch
    Nx   = size(bndx.value)[1] - 1
    Ny   = size(bndy.value)[1] - 1
    bval = zeros(Nx+1, Ny+1)
    if bndx.value[1] != bndy.value[1]
        error("Inconsistent boundary conditions.")
    else
        bval[:, 1] = bndx.value
        bval[1, :] = bndy.value 
    end
    return Patch(loc, shapeL2H(shapeH2L(derivOP + boundaryOP) \ shapeH2L(RHS + bval), Nx, Ny))  
end

function extractPatchCoeffs(patch::Patch)::Array{Float64,2}
    fnodal  = patch.value
    Nx      = size(fnodal)[1] - 1
    Ny      = size(fnodal)[2] - 1
    invndmx = inv(vandermonde(Nx))
    invndmy = inv(vandermonde(Ny))
    fmodal  = zeros(Nx+1, Ny+1)
    for index in CartesianRange(size(fmodal))
        i = index.I[1]
        j = index.I[2]
        fmodal[i,j] = sum(invndmx[i,m]*fnodal[m,n]*invndmy[j,n] for m in 1:Nx+1, n in 1:Ny+1) 
    end
    return fmodal
end

function interpolatePatch(patch::Patch, Nx::Int, Ny::Int)::Patch
    fmodal = extractPatchCoeffs(patch)
    Nxp    = size(fmodal)[1] - 1
    Nyp    = size(fmodal)[2] - 1
    pvndmx = pseudovandermonde(Nxp, chebgrid(Nx))
    pvndmy = pseudovandermonde(Nyp, chebgrid(Ny))
    fnodal = zeros(Nx+1, Ny+1)
    for index in CartesianRange(size(fnodal))
        i = index.I[1]
        j = index.I[2]
        fnodal[i,j] = sum(pvndmx[i,m]*fmodal[m,n]*pvndmy[j,n] for m in 1:Nxp+1, n in 1:Nyp+1)
    end
    return Patch(patch.loc, fnodal)
end
