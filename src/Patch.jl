#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

struct Patch
    loc::Array{Int,1}
    value::Array{Float64,2}
end

struct Boundary
    kind::Symbol
    value::Array{Float64,1}
end

function getPB(patch::Patch, s::Symbol)::Boundary
    if s==:R
        boundary = Boundary(:R, patch.value[end, :])
    elseif s==:C
        boundary = Boundary(:C, patch.value[:, end])
    else 
        error("Unknown symbol passed")
    end
    return boundary
end

function calcPatch(loc::Array{Int,1}, bnd0::Boundary, bnd1::Boundary, operator::Array{Float64, 4})::Patch
    N = size(operator)[1] - 1
    B = zeros(N+1, N+1)
    if bnd0.value[1] != bnd1.value[1]
        error("Inconsistent boundary conditions.")
    else
        B[1, :] = bnd0.value
        B[:, 1] = bnd1.value
    end
    return Patch(loc, shapeB(reshapeA(operator) \ reshapeB(B))) 
end

# XXX: Need to generalize these for more than one dimension. Implement index-wise multiplication
function prolongation2D(patch::Patch, M::Int)::Dict{Array{Int,1}, Patch}
    # Go from N modes in the entire patch to N modes in each subpatch
    N  = size(patch.value)[1] - 1
    xg = Float64[chebx(i,N) for i in 1:N+1]
    vx = vy = vandermonde(N, xg)
    dbase = Dict()
    for k in 1:M, l in 1:M
        loc = [k,l] 
        xp = Float64[coordtrans(M, [chebx(i,N),chebx(1,N)], loc)[1] for i in 1:N+1] 
        yp = Float64[coordtrans(M, [chebx(1,N),chebx(j,N)], loc)[2] for j in 1:N+1]
        px = vandermonde(N, xp) 
        py = vandermonde(N, yp)
        fgrid  = patch.value
        fpatch = px*inv(vx)*fgrid*inv(vy')*py'
        dbase[loc] = Patch(loc, fpatch) 
    end
    return dbase
end

# XXX: Need to generalize this to more dimensions.
# Can you do the restriction patch-wise?
function restriction2D(dbase::Dict{Array{Int, 1}, Patch}, M::Int)::Patch
    # Go from N modes in each patch to N modes in the entire patch
    N = size(dbase[[1,1]].value)[1] - 1
    fcpatches = zeros((N+1)*M, (N+1)*M)
    cPx = zeros((N+1)*M, (N+1)*M)
    cPy = zeros((N+1)*M, (N+1)*M)
    xg  = Float64[chebx(i,N) for i in 1:N+1]
    vx  = vy = vandermonde(N, xg)
    for i in 1:M, j in 1:M
        li = 1+(i-1)*(N+1)
        lj = 1+(j-1)*(N+1)
        loc = [i,j]
        xp = Float64[coordtrans(M, [chebx(i,N),chebx(1,N)], loc)[1] for i in 1:N+1]
        yp = Float64[coordtrans(M, [chebx(1,N),chebx(j,N)], loc)[2] for j in 1:N+1]
        px = vandermonde(N, xp)
        py = vandermonde(N, yp)
        cPx[li:li+N, lj:lj+N]   = px*inv(vx) 
        cPy[li:li+N, lj:lj+N]   = inv(vy')*py'
        fcpatches[li:li+N, lj:lj+N] = dbase[loc].value   
    end
    fglobalgrid = pinv(cPx)*fcpatches*pinv(cPy)
    return Patch([1,1], fglobalgrid) 
end

function prolongation1D(fxgrid::Array{Float64,1}, M::Int)::Dict{Int, Array{Float64,1}}
    # Function on a single patch >  dictionary with function vals on smaller subpatches
    # Goes from N modes in the global patch to N modes in each of the individual patches
    N     = size(fxgrid)[1] - 1
    xg    = Float64[chebx(i, N) for i in 1:N+1]
    vndm  = vandermonde(N,xg)
    dbase = Dict()
    for i in 1:M
        loc  = [1, i]
        xp   = Float64[coordtrans(M, [chebx(i,N),chebx(1,N)], loc)[1] for i in 1:N+1]
        px   = vandermonde(N,xp)
        P    = px*inv(vndm)
        dbase[i] = P*fxgrid
    end 
    return dbase
end     

function restriction1D(dbase::Dict{Int, Array{Float64,1}}, M::Int)::Array{Float64,1}
    # Function on multiple subpatches (dictionary) > Function on a single patch 
    # Go from N modes in each subpatch to  N modes in one single whole patch
    N    = size(dbase[1])[1] - 1
    xg   = Float64[chebx(i, N) for i in 1:N+1]
    vndm = vandermonde(N,xg)
    P    = Float64[]
    fxp  = Float64[]
    # loop over each patch, joining the P's and the datasets 
    for i in 1:M
        loc  = [1,i]
        xp   = Float64[coordtrans(M, [chebx(i,N),chebx(1,N)], loc)[1] for i in 1:N+1]
        px   = vandermonde(N,xp)
        p2P  = px*inv(vndm)
        fxp  = vcat(fxp, dbase[i])
        P    = vcat(P, p2P)
    end 
    # take the pseudo-inverse of the prolongation operator
    fxglobal = pinv(P)*fxp
    return fxglobal
end 

