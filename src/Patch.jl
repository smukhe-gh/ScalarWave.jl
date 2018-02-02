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

function pushforward(patch::Patch, loc::Array{Int,1})::Patch
    N  = size(patch.value) - 1
    x  = Float64[chebx(i,N) for i in 1:N+1]
    vx = vy = vandermonde(N, x)
    xp = Float64[coordtrans(M, [chebx(i,N),chebx(1,N)], loc)[1] for i in 1:N+1] 
    yp = Float64[coordtrans(M, [chebx(1,N),chebx(j,N)], loc)[2] for j in 1:N+1]
    px = vandermonde(N, xp) 
    py = vandermonde(N, xp)
    
    fgrid  = patch.value
    fpatch = px*inv(vx)*fgrid*inv(vy')*py'
    return Patch(loc, fpatch)
end

function pullback(N::Int, M::, dbase::Dict{Array{Int, 1}, Patch})::Patch
    # XXX: Can you do the restriction patch-wise?
    NG    = (N+1)*M
    fgrid = zeros(NG, NG)
    x  = Float64[chebx(i,NG) for i in 1:NG+1]
    vx = vy = vandermonde(NG, x)

    for i in 1:M, j in 1:M
        li = 1+(i-1)*(N+1)
        lj = 1+(j-1)*(N+1)
        xp = Float64[coordtrans(M, [chebx(i,N),chebx(1,N)], loc)[1] for i in 1:N+1]
        yp = Float64[coordtrans(M, [chebx(1,N),chebx(j,N)], loc)[2] for j in 1:N+1]
        px = vandermonde(N, xp)
        py = vandermonde(N, xp)
        fpatch = dbase[[i,j]]
        fgrid[li:li+N, lj:lj+N] = pinv(px*inv(vx))*fpatch*pinv(inv(vy')*py')   
    end
    return Patch([0,0], fgrid) 
end
