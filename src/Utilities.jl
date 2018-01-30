#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function delta{T<:Int}(i::T, j::T)::Float64
	return i==j ? 1 : 0
end

function coordtrans{T<:Int}(M::T, point::Array{Float64, 1}, loc::Array{T, 1})::Array{Float64, 1}
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
    return Float64[cheb(m,x) for m in 0:N, x in nodes]
end

function pconvergence(N::Int)::Float64
    dbase = distribute(N, 1, x-> sin(pi*x), y->sin(pi*y))
    chebGridData    = dbase[[1,1]]
    gaussGrid, w    = gauss(2*N)
    gaussGrid       = -gaussGrid   # flipping the nodes to be consistent with chebgrid
    exactGridData   = Float64[sin(pi*i) + sin(pi*j) for i in gaussGrid, j in gaussGrid]
    interpGridData  = interpolatePatch(chebGridData, gaussGrid, gaussGrid).value
    errorGridData   = interpGridData - exactGridData
    L2errorGridData = sqrt((w'*(errorGridData.^2)*w)/(w'*(exactGridData.^2)*w))
    return L2errorGridData
end

function hconvergence(M::Int)::Float64
    dbase = distribute(12, M, x-> sin(pi*x), y->sin(pi*y))
    gaussGrid, w  = gauss(12*M)
    gaussGrid     = -gaussGrid     # flipping the nodes to be consistent with chebgrid
    gaussGridData = Float64[sin(pi*i) + sin(pi*j) for i in gaussGrid, j in gaussGrid]
    chebGridData  = zeros(12*M, 12*M)
    for i in 1:M, j in 1:M
        li = 1+(i-1)*12
        lj = 1+(j-1)*12
        gaussLocalGridy = gaussGrid[li:li+11]
        gaussLocalGridx = gaussGrid[lj:lj+11]
        chebPatchData   = dbase[[i,j]]
        interpPatchData  = interpolatePatch(chebPatchData, gaussLocalGridx, gaussLocalGridy).value
        chebGridData[li:li+11, lj:lj+11] = interpPatchData
    end
    errorGridData   = chebGridData - gaussGridData
    L2errorGridData = sqrt((w'*(errorGridData.^2)*w)/(w'*(gaussGridData.^2)*w))
    return L2errorGridData
 end

