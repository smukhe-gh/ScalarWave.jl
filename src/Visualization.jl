 #--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 02-2018
#--------------------------------------------------------------------

function rect(x,y,w,h)
    return Shape(x + [0,w,w,0], y + [0,0,h,h]) 
end

function plotpatch(patch::Patch, M::Int)
    N  = size(patch)[1] - 1
    xp = Float64[coordtrans(M, [chebx(i,N),chebx(1,N)], loc)[1] for i in 1:N+1]
    yp = Float64[coordtrans(M, [chebx(i,N),chebx(1,N)], loc)[1] for i in 1:N+1]
    wx = Float64[coordtrans(M, [wx, 0], loc)[1] for wx in cumsum(chebweights(N))]
    wy = Float64[coordtrans(M, [0, wy], loc)[2] for wy in cumsum(chebweights(N))]
end
