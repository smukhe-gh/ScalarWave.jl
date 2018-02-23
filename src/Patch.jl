#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function getPatchIC{T<:Integer}(Nx::T, Ny::T, M::T, loc::Array{T,1}, fn::Function, s::Int)::Boundary
    if s == 0
        xg = chebgrid(Nx, M, loc[1])
        return Boundary(0, fn.(xg))
    elseif s == 1
        yg = chebgrid(Ny, M, loc[2])
        return Boundary(1, fn.(yg))
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

function calcPatch(loc::Array{Int,1}, bndx::Boundary, bndy::Boundary, operator::Array{Float64, 4})::Patch
    Nx   = size(operator)[1] - 1
    Ny   = size(operator)[2] - 1
    bval = zeros(Nx+1, Ny+1)
    if bnd0.value[1] != bnd1.value[1]
        error("Inconsistent boundary conditions.")
    else
        bval[:, 1] = bndx.value
        bval[1, :] = bndy.value 
    end
    return Patch(loc, shapeL2H(shapeH2L(operator) \ shapeH2L(bval))) 
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

function projectonPatchBnd(fn::Function, Nx::Int)::Array{Float64,1}
    coeffs = zeros(Nx+1)
    fonbnd = zeros(Nx+1)
    vndmx  = vandermonde(Nx)
    for m in 0:Nx
        integ = quadgk(x->fn(cos(x))*cos(m*x), 0, pi)[1]
        (m == 0) ? coeffs[m+1] = integ/pi : coeffs[m+1] = integ/(pi/2)
    end
    for nx in 1:Nx+1
        fonbnd[nx] = sum(vndmx[nx,j]*coeffs[j] for j in 1:Nx+1)
    end
    return fonbnd 
end

function projectonPatchBnd(fn::Function, Nx::Int, M::Int, loc::Int)::Array{Float64,1}
    coeffs = zeros(Nx+1)
    fonbnd = zeros(Nx+1)
    vndmx  = vandermonde(Nx)
    for m in 0:Nx
        integ = quadgk(xp->fn(cos(coortransL2G(M, loc, xp)))*cos(m*xp), 0, pi)[1]
        (m == 0) ? coeffs[m+1] = integ/pi : coeffs[m+1] = integ/(pi/2)
    end
    for nx in 1:Nx+1
        fonbnd[nx] = sum(vndmx[nx,j]*coeffs[j] for j in 1:Nx+1)
    end
    return fonbnd 
end

function projectonPatch(fn::Function, Nx::Int, Ny::Int, M::Int, loc::Array{Float64,1})::Array{Float64,2}
    fmodal   = zeros(Nx+1, Ny+1)
    fonpatch = zeros(Nx+1, Ny+1)
    vndmx    = vandermonde(Nx)
    vndmy    = vandermonde(Ny)
    for m in 0:Nx,  n in 0:Ny
        integxy = hcubature(xp->fn(cos(coordtransL2G(M,loc[1],xp[1])), cos(coordtrans(M,loc[2],xp[2])))*cos(m*x[1])*cos(n*x[2]), (0, 0), (pi, pi))
        if m == n == 0
            fmodal[m,n] = integxy/(pi^2)
        elseif m == 0 || n == 0
            fmodal[m,n] = integxy/(pi^2/2)
        else
            fmodal[m,n] = integxy/(pi^2/4)
        end
    end
    for index in CartesianRange(size(fonpatch))
         i = index.I[1]
         j = index.I[2]
         fonpatch[i,j] = sum(vndmx[i,m]*fmodal[m,n]*vndmy[j,n] for m in 1:Nx+1, n in 1:Ny+1)
    end 
    return fonpatch
end
