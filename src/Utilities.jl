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

function L2norm{T<:Array{Float64,2}}(errorGridData::T, exactGridData::T, w::Array{Float64,1})::Float64
    return sqrt((w'*(errorGridData.^2)*w)/(w'*(exactGridData.^2)*w))
end

function L1norm{T<:Array{Float64,2}}(errorGridData::T)::Float64
    return maximum(abs.(errorGridData))
end

function chebgrid(N::Int)::Array{Float64,1}
    return Float64[chebx(i,N) for i in 1:N+1] 
end

function pconvergence(N::Int, M::Int)::Float64
    fdbase  = distribute(N, 1, x-> sin(pi*x), y->sin(pi*y))
    exactGrid = Float64[sin(pi*i) + sin(pi*j) for i in chebgrid(N), j in chebgrid(N)]
    # project both the computed and the exact solution onto 4 patches
    fcompute  = prolongation2D(fdbase[[1,1]], M)
    fexact    = prolongation2D(Patch([1,1], exactGrid), M)
    errornorm = zeros(M*M)
    # compute the error patch-wise
    for i in 1:M, j in 1:M
        errornorm[i+j] = L1norm(fexact[[i,j]].value - fcompute[[i,j]].value)
    end
    return maximum(errornorm)
end

function hconvergence(N::Int, M::Int)::Float64
    fcompute  = distribute(N, M, x-> sin(pi*x), y->sin(pi*y))
    exactGrid = Float64[sin(pi*i) + sin(pi*j) for i in chebgrid(N), j in chebgrid(N)]
    # project the exact solution on the same number of patches as the computation
    fexact    = prolongation2D(Patch([1,1], exactGrid), M)
    errornorm = zeros(M*M)
    # compute the error patch-wise
    for i in 1:M, j in 1:M
        errornorm[i+j] = L1norm(fexact[[j,i]].value - fcompute[[i,j]].value)
    end
    return maximum(errornorm)
end
