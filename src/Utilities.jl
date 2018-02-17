#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function delta{T<:Int}(i::T, j::T)::Float64
	return i==j ? 1 : 0
end

function coordtrans{T<:Int}(M::T, point::Array{Float64, 1}, loc::Array{T, 1})::Array{Float64, 1}
    if maximum(loc) > M
        error("Location incompatible with the number of patches")
    end
    x, y = point
    s = Float64[d for d in M-1:-2:-M+1]
	xp = (x + s[loc[2]])/M
	yp = (y + s[loc[1]])/M
	return [xp,yp]
end

function reshapeA(op4N::Array{Float64,4})::Array{Float64,2}
    N = size(op4N)[1] - 1
    op2N = reshape(op4N, ((N+1)^2, (N+1)^2))
    return op2N
end

function shapeA(op2N::Array{Float64,2})::Array{Float64,4}
    N = Int(sqrt(size(op2N)[1])) - 1
    op4N = reshape(op2N, (N+1, N+1, N+1, N+1))
    return op4N
end

function reshapeB(b2N::Array{Float64,2})::Array{Float64,1}
    N = size(b2N)[1] - 1
    b1N = reshape(b2N, (N+1)^2)
    return b1N
end

function shapeB(b1N::Array{Float64,1})::Array{Float64,2}
    N = Int(sqrt(size(b1N)[1])) - 1
    b2N = reshape(b1N, (N+1, N+1))
    return b2N
end

function vandermonde(N::Int, nodes::Array{Float64, 1})::Array{Float64,2}
    return Float64[cheb(m,x) for x in nodes, m in 0:N]
end

# XXX: Test cases missing from here. 

function L2norm{T<:Array{Float64,2}}(errorGridData::T, exactGridData::T, w::Array{Float64,1})::Float64
    return sqrt((w'*(errorGridData.^2)*w)/(w'*(exactGridData.^2)*w))
end

function LInfnorm{T<:Array{Float64,2}}(errorGridData::T)::Float64
    return maximum(abs.(errorGridData))
end

function chebgrid(N::Int)::Array{Float64,1}
    return Float64[chebx(i,N) for i in 1:N+1] 
end

function chebweights(N::Int)::Array{Float64,1}
    return Float64[chebw(i,N) for i in 1:N+1]
end

# TODO: One can test the convergence routines separately, just by representing functions on grids.
function pconvergence(N::Int, fbnd1:: Function, fbnd2::Function, fsol::Function)::Float64
    dbase = distribute(N, 1, fbnd1, fbnd2)
    chebGrid        = Float64[chebx(i, N) for i in 1:N+1]
    chebGridData    = dbase[[1,1]]
    gaussGrid, w    = gauss(2*N)
    gaussGrid       = - gaussGrid
    exactGridData   = Float64[fsol(i,j) for i in gaussGrid, j in gaussGrid]
    interpGridData  = interpolatePatch(chebGridData, gaussGrid, gaussGrid).value
    errorGridData   = interpGridData - exactGridData
    L2errorGridData = sqrt((w'*(errorGridData.^2)*w)/(w'*(exactGridData.^2)*w))
    return L2errorGridData
end

# TODO: One can test the convergence routines separately, just by representing functions on grids.
function hconvergence(N::Int, M::Int, fbnd1:: Function, fbnd2::Function, fsol::Function)::Float64
    dbase = distribute(N, M, fbnd1, fbnd2)
    chebGrid        = Float64[chebx(i, N) for i in 1:N+1]
    gaussGrid, w    = gauss(2*N)
    gaussGrid       = -gaussGrid
    errorvec        = zeros(M*M)
    for i in 1:M, j in 1:M
        loc             = [i,j]
        gaussLocalx     = Float64[coordtrans(M, [x, 0], loc)[1] for x in gaussGrid]
        gaussLocaly     = Float64[coordtrans(M, [0, y], loc)[2] for y in gaussGrid] 
        exactGridData   = Float64[fsol(x,y) for x in gaussLocalx, y in gaussLocaly]
        interpGridData  = interpolatePatch(dbase[loc], gaussGrid, gaussGrid).value
        errorGridData   = interpGridData - exactGridData
        L2errorGridData = sqrt((w'*(errorGridData.^2)*w)/(w'*(exactGridData.^2)*w))
        errorvec[i+j]   = L2errorGridData
    end
    return sum(errorvec)
end

# TODO: Add the option to store params.
function savegrid(dbase::Dict{Array{Int,1}, Patch})
    datetime =  DateTime(now())
    jldopen("$datetime.jld", "w") do file
        addrequire(file, ScalarWave)
        write(file, "patches", dbase)
    end
end

function projectboundary(func::Function, N::Int)::Array{Float64,1}
    coeffs = zeros(N+1)
    for m in 0:N
        if m == 0
            # TODO: Test guassian integration routines
            coeffs = quadgk(x->func(cos(x))*cos(m*x), 0, pi; abstol=0, maxevals=10^7, order=2*N, norm=vecnorm)/pi
        else
            coeffs = quadgk(x->func(cos(x))*cos(m*x), 0, pi; abstol=0, maxevals=10^7, order=2*N, norm=vecnorm)/(pi/2)
        end
    end
    # TODO: Make this operation work element-wise. Can you put this inside the upper loop?
    return vandermonde(N, chebgrid(N))*coeffs
end
